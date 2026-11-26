#!/usr/bin/env python3
"""Add county-specific Census PLACE populations to Ohio's existing five-city graph.

Reverse-engineer (and verify) the five-city contraction using MEMBER_NAMES, then
produce new, non-overwriting copies of both graphs with PLACE_COUNTY dictionaries.
Preserve all original nodes, edges, weights, coordinates, labels and metadata.

Run from the Ohio repository:

  python3 build_oh_place_county_graph.py

If the July original graph uses PLACE names but the November graph uses GEOIDs,
provide the Census PLACE polygons as a crosswalk source:

  python3 build_oh_place_county_graph.py \
      --census-places data/raw/tl_2024_39_place

Alternatively provide a CSV with NAME and GEOID columns:

  python3 build_oh_place_county_graph.py --place-crosswalk-csv name_geoid.csv

Only standard-library Python is needed when keys already match. GeoPandas is
only needed to read a Census polygon file when the keys differ.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
import os
from pathlib import Path
import re
import tempfile

CITIES = frozenset(("Akron", "Cincinnati", "Cleveland", "Dayton", "Toledo"))


def require(condition: bool, detail: str) -> None:
    if not condition:
        raise ValueError(detail)


def integer(value, description: str) -> int:
    require(not isinstance(value, bool), f"Invalid integer for {description}: {value!r}")
    try:
        n = int(value)
        valid = n == float(value)
    except (TypeError, ValueError, OverflowError):
        valid = False
        n = 0
    require(valid and n >= 0, f"Invalid nonnegative integer for {description}: {value!r}")
    return n


def check_place_dict(value, label: str) -> dict[str, int]:
    require(isinstance(value, dict), f"{label} must be a dictionary, got {type(value).__name__}")
    out = {}
    for key, val in value.items():
        require(isinstance(key, str) and key.strip() and "|" not in key,
                f"Invalid PLACE identifier {key!r} at {label}")
        pop = integer(val, f"{label}[{key!r}]")
        if pop:
            out[key] = pop
    return dict(sorted(out.items()))


def node_index(graph: dict, label: str) -> tuple[dict, dict]:
    require(isinstance(graph, dict) and isinstance(graph.get("nodes"), list),
            f"{label}: expected JSON graph object with a nodes list")
    require(graph.get("directed", False) is False, f"{label}: directed graph not supported")
    require(graph.get("multigraph", False) is False, f"{label}: multigraph not supported")
    by_id = {}
    by_name = {}
    for node in graph["nodes"]:
        require(isinstance(node, dict) and "id" in node and "NAME" in node,
                f"{label}: node missing id or NAME")
        nid, name = node["id"], node["NAME"]
        require(isinstance(nid, (int, str)) and not isinstance(nid, bool),
                f"{label}: invalid node id {nid!r}")
        require(isinstance(name, str) and name, f"{label}: invalid NAME {name!r}")
        require(nid not in by_id, f"{label}: duplicate id {nid!r}")
        require(name not in by_name, f"{label}: duplicate NAME {name!r}")
        by_id[nid] = node
        by_name[name] = node
    return by_id, by_name


def edges_of(graph: dict, by_id: dict, label: str) -> set[frozenset]:
    """Read undirected edge sets in NetworkX node-link or adjacency JSON formats."""
    edges: set[frozenset] = set()
    if "adjacency" in graph:
        adjacency = graph["adjacency"]
        require(isinstance(adjacency, list) and len(adjacency) == len(graph["nodes"]),
                f"{label}: adjacency length differs from nodes")
        for node, neighbors in zip(graph["nodes"], adjacency):
            require(isinstance(neighbors, list), f"{label}: invalid adjacency entry")
            source = node["id"]
            for neighbor in neighbors:
                require(isinstance(neighbor, dict) and "id" in neighbor,
                        f"{label}: expected adjacency neighbors with id")
                target = neighbor["id"]
                require(target in by_id, f"{label}: edge target {target!r} not in nodes")
                require(source != target, f"{label}: original graph contains a self-loop at {source!r}")
                edges.add(frozenset((source, target)))
    else:
        key = "links" if "links" in graph else "edges" if "edges" in graph else None
        require(key is not None, f"{label}: cannot locate graph edges (adjacency/links/edges)")
        require(isinstance(graph[key], list), f"{label}: {key} is not a list")
        for item in graph[key]:
            require(isinstance(item, dict) and "source" in item and "target" in item,
                    f"{label}: {key} item has no source/target")
            source, target = item["source"], item["target"]
            require(source in by_id and target in by_id,
                    f"{label}: unknown edge endpoint {source!r},{target!r}")
            require(source != target, f"{label}: self-loop at {source!r}")
            edges.add(frozenset((source, target)))
    require(edges, f"{label}: no edges found")
    return edges


def crosswalk_from_csv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames is not None, f"{path}: empty CSV")
        require({"NAME", "GEOID"}.issubset(reader.fieldnames),
                f"{path}: CSV needs NAME and GEOID headers")
        pairs = [(row["NAME"].strip(), row["GEOID"].strip()) for row in reader]
    return validate_crosswalk(pairs, path)


def crosswalk_from_census(path: Path) -> dict[str, str]:
    try:
        import geopandas as gpd
    except ImportError as exc:
        raise RuntimeError("GeoPandas required for --census-places; "
                           "use --place-crosswalk-csv NAME,GEOID instead") from exc
    gdf = gpd.read_file(path)
    require({"NAME", "GEOID"}.issubset(gdf.columns),
            f"{path}: Census data must contain NAME and GEOID fields")
    return validate_crosswalk(zip(gdf["NAME"], gdf["GEOID"]), path)


def validate_crosswalk(pairs, source: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for name, geoid in pairs:
        name, geoid = str(name).strip(), str(geoid).strip()
        require(name and re.fullmatch(r"39\d{5}", geoid) is not None,
                f"{source}: invalid Census PLACE NAME/GEOID: {name!r},{geoid!r}")
        if name in mapping:
            require(mapping[name] == geoid, f"{source}: ambiguous PLACE name {name!r}")
        mapping[name] = geoid
    require(mapping, f"{source}: no PLACE crosswalk records")
    return mapping


def reconcile_places(orig_by_name: dict, sampled_by_name: dict,
                     crosswalk: dict[str, str] | None) -> tuple[dict[str, dict[str, int]], str]:
    """Normalize original PLACE entries to the exact identifiers in sampling JSON."""
    raw = {name: check_place_dict(node.get("PLACE"), f"original PLACE of {name}")
           for name, node in orig_by_name.items()}
    converted = {}
    for name, place in raw.items():
        new = defaultdict(int)
        for key, pop in place.items():
            if crosswalk and not re.fullmatch(r"39\d{5}", key):
                require(key in crosswalk, f"Original PLACE {key!r} lacks Census crosswalk entry")
                key = crosswalk[key]
            new[key] += pop
        converted[name] = dict(sorted(new.items()))
    # Compare all uncontracted nodes. If crosswalk missing, signal actionable error.
    mismatches = []
    for name, node in sampled_by_name.items():
        if "PROTECTED_CITY" not in node:
            sampled_place = check_place_dict(node.get("PLACE"), f"sampling PLACE of {name}")
            if converted.get(name) != sampled_place:
                mismatches.append((name, converted.get(name), sampled_place))
    if mismatches:
        name, original, sampling = mismatches[0]
        original_keys = sorted(k for places in raw.values() for k in places)[:8]
        sampled_keys = sorted(k for node in sampled_by_name.values()
                              for k in check_place_dict(node.get("PLACE"),
                                                        f"sampling PLACE of {node['NAME']}"))[:8]
        raise ValueError(
            f"PLACE representation/population mismatch in {len(mismatches)} uncontracted nodes.\n"
            f"First mismatch {name!r}: original={original!r}; sampling={sampling!r}.\n"
            f"Original key examples: {original_keys!r}; sampling key examples: {sampled_keys!r}.\n"
            "If the original keys are NAMEs but sampling uses GEOIDs, rerun with "
            "--census-places data/raw/tl_2024_39_place or --place-crosswalk-csv. "
            "If already using a crosswalk, STOP: underlying PLACE populations differ."
        )
    return converted, "crosswalk" if crosswalk else "original_keys"


def make_city_member_mapping(original_by_name: dict, sampled: dict):
    city_nodes = {}
    mapping_by_original = {}
    city_reports = {}
    for node in sampled["nodes"]:
        city = node.get("PROTECTED_CITY")
        if city is not None:
            require(city in CITIES, f"Unexpected PROTECTED_CITY {city!r}")
            require(city not in city_nodes, f"Repeated protected-city node for {city}")
            members = node.get("MEMBER_NAMES")
            require(isinstance(members, list) and members and all(isinstance(x, str) for x in members),
                    f"City {city}: MEMBER_NAMES should be a nonempty list of precinct names")
            require(len(members) == len(set(members)), f"City {city}: duplicate MEMBER_NAMES")
            city_nodes[city] = node
        else:
            require("MEMBER_NAMES" not in node,
                    f"Noncity node {node['NAME']!r} unexpectedly has MEMBER_NAMES")
            members = [node["NAME"]]
        for name in members:
            require(name in original_by_name, f"Sampling node {node['NAME']!r} refers to unknown precinct {name!r}")
            require(name not in mapping_by_original,
                    f"Original precinct {name!r} appears in multiple contracted nodes")
            mapping_by_original[name] = node
        if city is not None:
            county_pop = Counter()
            for name in members:
                original = original_by_name[name]
                county_pop[original["COUNTY"]] += integer(original["POP20"], f"POP20 in {name}")
            city_reports[city] = {
                "node": node["NAME"],
                "member_count": len(members),
                "precinct_pop20_by_county": dict(sorted(county_pop.items())),
                "pop20": sum(county_pop.values()),
            }
    require(set(city_nodes) == CITIES,
            f"Expected precisely {sorted(CITIES)}, found {sorted(city_nodes)}")
    require(set(mapping_by_original) == set(original_by_name),
            f"Missing {len(set(original_by_name)-set(mapping_by_original))} original precincts")
    return mapping_by_original, city_reports


def topology_check(original: dict, sampled: dict, orig_by_id: dict,
                   samp_by_id: dict, mapping_by_original: dict) -> dict:
    original_edges = edges_of(original, orig_by_id, "original")
    sampled_edges = edges_of(sampled, samp_by_id, "five-city")
    expected = set()
    collapsed_internal = 0
    for edge in original_edges:
        first, second = tuple(edge)
        a = mapping_by_original[orig_by_id[first]["NAME"]]["id"]
        b = mapping_by_original[orig_by_id[second]["NAME"]]["id"]
        if a == b:
            collapsed_internal += 1
        else:
            expected.add(frozenset((a, b)))
    if expected != sampled_edges:
        missing = expected - sampled_edges
        added = sampled_edges - expected
        def example(values):
            return [sorted(map(str, edge)) for edge in list(values)[:5]]
        raise ValueError(
            f"Five-city topology does NOT equal the graph obtained by contracting original edges. "
            f"Missing {len(missing)} edges, extra {len(added)} edges. "
            f"Missing examples={example(missing)}, extra examples={example(added)}. "
            "STOP; the saved graph was not made by this simple contraction alone."
        )
    return {
        "original_edges": len(original_edges),
        "sampling_edges": len(sampled_edges),
        "internal_edges_removed": collapsed_internal,
        "topology_exact": True,
    }


def sum_dicts(items):
    totals = Counter()
    for mapping in items:
        totals.update(mapping)
    return dict(sorted((key, val) for key, val in totals.items() if val))


def enrich(original: dict, sampled: dict, crosswalk: dict | None,
           expected_original: int = 7928, expected_sampled: int = 6975,
           expected_counties: int = 88) -> dict:
    orig_by_id, orig_by_name = node_index(original, "original")
    samp_by_id, samp_by_name = node_index(sampled, "five-city")
    require(len(orig_by_name) == expected_original,
            f"Expected {expected_original} original nodes, got {len(orig_by_name)}")
    require(len(samp_by_name) == expected_sampled,
            f"Expected {expected_sampled} sampling nodes, got {len(samp_by_name)}")

    mapping_by_original, city_reports = make_city_member_mapping(orig_by_name, sampled)
    original_place, key_mode = reconcile_places(orig_by_name, samp_by_name, crosswalk)
    topology = topology_check(original, sampled, orig_by_id, samp_by_id, mapping_by_original)

    # With an explicit Census name->GEOID crosswalk, independently establish
    # that every precinct containing each protected city's population belongs
    # to that city's protected supernode. This is stricter than simply
    # retaining the five PREEXISTING contracted nodes.
    city_membership_audit = {"checked": crosswalk is not None}
    if crosswalk is not None:
        for city in sorted(CITIES):
            require(city in crosswalk, f"Census crosswalk missing protected city {city!r}")
            geoid = crosswalk[city]
            positive = [name for name, p in original_place.items() if p.get(geoid, 0) > 0]
            require(positive, f"No original precinct population identified for {city} ({geoid})")
            wrong = [name for name in positive
                     if mapping_by_original[name].get("PROTECTED_CITY") != city]
            require(not wrong,
                    f"The saved contraction does not preserve all {city} population: "
                    f"{len(wrong)} original precincts carrying {geoid} are outside "
                    f"its protected city node; first={wrong[:8]}")
            city_membership_audit[city] = {
                "geoid": geoid,
                "original_precincts_with_positive_place_population": len(positive),
                "whole_city_in_protected_node": True,
            }

    county_names = {node["COUNTY"] for node in original["nodes"]}
    require(len(county_names) == expected_counties,
            f"Expected {expected_counties} counties, got {len(county_names)}")
    require(all(isinstance(c, str) and c and "|" not in c for c in county_names),
            "COUNTY values must be nonempty strings without '|'")

    original_county_place = {}
    total_place = 0
    for name, node in orig_by_name.items():
        pop = integer(node["POP20"], f"original POP20 of {name}")
        place = original_place[name]
        require(sum(place.values()) <= pop,
                f"Original precinct {name!r} has PLACE population > POP20")
        county = node["COUNTY"]
        county_place = {f"{pid}|{county}": p for pid, p in place.items()}
        require("PLACE_COUNTY" not in node,
                f"Original graph already has PLACE_COUNTY on {name}")
        node["PLACE_COUNTY"] = county_place
        original_county_place[name] = county_place
        total_place += sum(place.values())

    city_member_sets = {city: [] for city in CITIES}
    total_sampling_place = 0
    total_sample_pop = 0
    for node in sampled["nodes"]:
        members = node["MEMBER_NAMES"] if "PROTECTED_CITY" in node else [node["NAME"]]
        member_pop = sum(integer(orig_by_name[name]["POP20"], f"original POP20 of {name}")
                         for name in members)
        sample_pop = integer(node["POP20"], f"sampling POP20 of {node['NAME']}")
        require(member_pop == sample_pop,
                f"{node['NAME']}: sampled POP20={sample_pop}, member sum={member_pop}")
        aggregated = sum_dicts(original_place[name] for name in members)
        sampled_place = check_place_dict(node.get("PLACE"),
                                         f"sampling PLACE of {node['NAME']}")
        require(aggregated == sampled_place,
                f"Sampling PLACE differs from sum of constituent precincts at "
                f"{node['NAME']!r}: original members={aggregated}, sampled={sampled_place}. "
                "STOP; may require different PLACE source data.")
        new_attr = sum_dicts(original_county_place[name] for name in members)
        # Strong local conservation: dropping county from the composite key
        # must recover every existing PLACE value, not merely its total.
        collapsed = Counter()
        for key, p in new_attr.items():
            place_id, _county = key.rsplit("|", 1)
            collapsed[place_id] += p
        require(dict(sorted(collapsed.items())) == sampled_place,
                f"PLACE_COUNTY does not reconstruct PLACE at {node['NAME']}")
        require("PLACE_COUNTY" not in node,
                f"Sampling graph already has PLACE_COUNTY on {node['NAME']}")
        node["PLACE_COUNTY"] = new_attr
        if "PROTECTED_CITY" in node:
            city_reports[node["PROTECTED_CITY"]]["place_county_populations"] = new_attr
        total_sample_pop += sample_pop
        total_sampling_place += sum(sampled_place.values())
        if "PROTECTED_CITY" in node:
            city_member_sets[node["PROTECTED_CITY"]] = members
    require(total_place == total_sampling_place,
            f"Statewide PLACE mismatch: original {total_place}, sampled {total_sampling_place}")
    original_pop = sum(integer(node["POP20"], f"original POP20 {node['NAME']}")
                       for node in original["nodes"])
    require(original_pop == total_sample_pop, "Statewide POP20 differs")
    county_units = sum_dicts(node["PLACE_COUNTY"] for node in sampled["nodes"])
    split_places = Counter()
    for key, population in county_units.items():
        place_id, _county = key.rsplit("|", 1)
        split_places[place_id] += 1
    for city, members in city_member_sets.items():
        require(members, f"Protected city {city} has no members")

    return {
        "original_nodes": len(orig_by_name),
        "five_city_nodes": len(samp_by_name),
        "protected_cities": city_reports,
        "protected_city_place_membership": city_membership_audit,
        "county_count": len(county_names),
        "place_key_normalization": key_mode,
        "statewide_pop20": original_pop,
        "statewide_place_population": total_place,
        "place_county_units": len(county_units),
        "places_in_multiple_counties": sum(count > 1 for count in split_places.values()),
        "topology": topology,
        "invariants": [
            "all original precincts represented exactly once",
            "when Census crosswalk is provided, all five CITY PLACE populations lie in their protected nodes",
            "exactly five protected city nodes retained",
            "sampled POP20 equals sum of original members for every node",
            "sampled PLACE equals sum of original members for every node",
            "PLACE_COUNTY sums to sampled PLACE for every node and PLACE",
            "all sampled graph edges equal quotient of original graph edges",
            "original node/link payloads left intact except new PLACE_COUNTY",
        ],
    }


def write_new_json(path: Path, obj) -> None:
    require(not path.exists(), f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = None
    try:
        with tempfile.NamedTemporaryFile("w", dir=path.parent,
                                         prefix=f".{path.name}.", suffix=".tmp",
                                         encoding="utf-8", delete=False) as out:
            temp = Path(out.name)
            json.dump(obj, out, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
            out.write("\n")
            out.flush()
            os.fsync(out.fileno())
        os.replace(temp, path)
    finally:
        if temp is not None:
            temp.unlink(missing_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--original", type=Path, default=Path("data/graph/OH_pct20.json"))
    ap.add_argument("--contracted", type=Path,
                    default=Path("data/graph/OH_pct20_preserved_five_cities.json"))
    ap.add_argument("--original-output", type=Path,
                    default=Path("data/graph/OH_pct20_place_county.json"))
    ap.add_argument("--contracted-output", type=Path,
                    default=Path("data/graph/OH_pct20_preserved_five_cities_place_county.json"))
    ap.add_argument("--report", type=Path,
                    default=Path("data/graph/OH_five_city_place_county_audit.json"))
    ap.add_argument("--census-places", type=Path, help="Census tl_2024_39_place shapefile/zip/path for NAME->GEOID")
    ap.add_argument("--place-crosswalk-csv", type=Path, help="CSV with NAME,GEOID columns")
    args = ap.parse_args()
    require(not (args.census_places and args.place_crosswalk_csv),
            "Use either --census-places or --place-crosswalk-csv, not both")
    inputs = [args.original.resolve(), args.contracted.resolve()]
    outputs = [args.original_output.resolve(), args.contracted_output.resolve(), args.report.resolve()]
    require(len(set(inputs + outputs)) == 5, "All input/output paths must be distinct")
    for path in inputs:
        require(path.is_file(), f"Missing input file: {path}")
    for path in outputs:
        require(not path.exists(), f"Refusing to overwrite existing file: {path}")
    with args.original.open(encoding="utf-8") as file:
        original = json.load(file)
    with args.contracted.open(encoding="utf-8") as file:
        sampled = json.load(file)
    if args.census_places:
        crosswalk = crosswalk_from_census(args.census_places)
    elif args.place_crosswalk_csv:
        crosswalk = crosswalk_from_csv(args.place_crosswalk_csv)
    else:
        crosswalk = None
    report = enrich(original, sampled, crosswalk)
    report["inputs"] = {
        "original": {"path": str(args.original)},
        "contracted": {"path": str(args.contracted)},
    }
    print(json.dumps(report, indent=2), flush=True)
    # Only write AFTER all membership, population, identity and topology checks pass.
    write_new_json(args.original_output, original)
    write_new_json(args.contracted_output, sampled)
    write_new_json(args.report, report)
    print("\nVALIDATED. Original input graphs were not modified.", flush=True)
    print("Original + PLACE_COUNTY:", args.original_output, flush=True)
    print("Five cities + PLACE_COUNTY:", args.contracted_output, flush=True)
    print("Audit:", args.report, flush=True)


if __name__ == "__main__":
    main()
