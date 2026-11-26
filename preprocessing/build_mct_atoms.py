#!/usr/bin/env python3
"""Construct block-population MCT atoms for the Ohio CycleWalk sampler.

No changes are made to the precinct graph. The block -> final graph-node
crosswalk is audited against *every* node's 2020 population before output.

Geography is kept in one 2024 TIGER vintage to match GeneratingJSON_OH.py.
Only LSAD 44 / CLASSFP T* townships and LSAD 25/47 / CLASSFP C*
incorporated cities and villages are selected. Thus CDPs and independent
municipal corporations in the COUSUB layer are NOT township units.

A municipality is nested only if its entire polygon (within a small geometry
sliver tolerance) lies inside one qualifying township. Other incorporated
places, including those spanning townships, belong to the cross group.

Output atoms aggregate 2020 census blocks by (graph node, township, place).
Township and place membership are computed independently for each block,
so rows form a DISJOINT population partition and admit joint diagnostics.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
import pandas as pd

OH_PLACE = re.compile(r"39\d{5}$")
OH_COUSUB = re.compile(r"39\d{8}$")
ATOM_FIELDS = ("node_name", "township_geoid", "place_geoid", "place_class", "population")


def census_id(value: object, pattern: re.Pattern[str], label: str) -> str:
    if pd.isna(value):
        raise ValueError(f"Missing {label} GEOID")
    key = str(value).strip()
    if not pattern.fullmatch(key):
        raise ValueError(f"Invalid {label} GEOID {key!r}; retain string leading zeroes")
    return key


def selected_geographies(source: str | Path, *, kind: str) -> gpd.GeoDataFrame:
    data = gpd.read_file(source)
    required = {"GEOID", "NAME", "LSAD", "CLASSFP", "geometry"}
    if not required.issubset(data.columns):
        raise ValueError(f"{kind} file missing columns {sorted(required - set(data.columns))}")
    if data.crs is None:
        raise ValueError(f"{kind} geometry has no CRS")
    if kind == "township":
        selected = data.loc[(data.LSAD.astype(str) == "44") &
                            data.CLASSFP.astype(str).str.startswith("T")].copy()
        pattern = OH_COUSUB
    elif kind == "place":
        selected = data.loc[data.LSAD.astype(str).isin(["25", "47"]) &
                            data.CLASSFP.astype(str).str.startswith("C")].copy()
        pattern = OH_PLACE
    else:
        raise ValueError(f"Unknown geographic type: {kind}")
    if selected.empty:
        raise ValueError(f"No qualifying {kind} units in {source}")
    selected["GEOID"] = [census_id(x, pattern, kind) for x in selected.GEOID]
    if selected.GEOID.duplicated().any():
        raise ValueError(f"Duplicate {kind} GEOIDs")
    if selected.geometry.isna().any() or selected.geometry.is_empty.any():
        raise ValueError(f"Empty {kind} geometries")
    if not selected.geometry.is_valid.all():
        raise ValueError(f"Invalid {kind} geometries; repair the source deliberately")
    return selected[["GEOID", "NAME", "CLASSFP", "LSAD", "geometry"]].reset_index(drop=True)


def place_nesting(places: gpd.GeoDataFrame, townships: gpd.GeoDataFrame,
                  *, relative_sliver: float = 1e-8,
                  absolute_sliver_m2: float = 1.0) -> pd.DataFrame:
    """Return fixed nested/cross labels independent of a districting plan.

    Uses true polygon difference, not representative-point membership.
    Each nested municipality must fit inside a SINGLE qualifying township.
    Ratios are diagnostic; small sliver tolerance only, never a 99% shortcut.
    """
    if relative_sliver < 0 or absolute_sliver_m2 < 0:
        raise ValueError("Sliver tolerances must be nonnegative")
    if places.crs != townships.crs:
        raise ValueError("PLACE and COUSUB layers must use identical CRS")
    if not places.crs.is_projected:
        raise ValueError("Area test requires projected coordinates")
    result = []
    tree = townships.sindex
    for p in places.itertuples(index=False):
        area = p.geometry.area
        if area <= 0:
            raise ValueError(f"Place {p.GEOID} has zero area")
        candidate_ids = tree.query(p.geometry, predicate="intersects")
        best_geoid, best_outside = "", area
        for j in candidate_ids:
            t = townships.iloc[int(j)]
            outside = p.geometry.difference(t.geometry).area
            if outside < best_outside:
                best_outside, best_geoid = outside, t.GEOID
        tolerance = max(absolute_sliver_m2, relative_sliver * area)
        is_nested = bool(best_geoid and best_outside <= tolerance)
        result.append({
            "place_geoid": p.GEOID,
            "place_name": p.NAME,
            "place_classfp": p.CLASSFP,
            "place_lsad": p.LSAD,
            "status": "nested" if is_nested else "cross",
            "parent_township_geoid": best_geoid if is_nested else "",
            "nearest_township_geoid": best_geoid,
            "outside_nearest_township_area_ratio": best_outside / area,
        })
    return pd.DataFrame(result).sort_values("place_geoid").reset_index(drop=True)


def graph_node_lookup(graph_nodes: list[dict]) -> tuple[dict, dict]:
    """Map underlying precinct components to the final post-articulation NAME.

    The production graph uses ORIG_NAME (possibly comma-joined), followed by
    __SD{senate}__HD{house}. Compare block-level enacted district labels so
    a precinct component cannot accidentally map to a different fragment.
    """
    mapping: dict[tuple[str, str, int, int], str] = {}
    graph_pop: dict[str, int] = {}
    for node in graph_nodes:
        name = str(node["NAME"])
        county = str(node["COUNTY"])
        sd = int(node["ENACTED_UPPER"])
        hd = int(node["ENACTED_LOWER"])
        pop = int(node["POP20"])
        if name in graph_pop or pop < 0:
            raise ValueError(f"Repeated graph name or negative population: {name}")
        graph_pop[name] = pop
        orig = str(node.get("ORIG_NAME", name.rsplit("__SD", 1)[0]))
        pieces = [x.strip() for x in orig.split(",")]
        if not all(pieces) or len(pieces) != len(set(pieces)):
            raise ValueError(f"Ambiguous ORIG_NAME on {name}: {orig}")
        for piece in pieces:
            key = (county, piece, sd, hd)
            if key in mapping:
                raise ValueError(f"Graph component appears in multiple nodes: {key}")
            mapping[key] = name
    return mapping, graph_pop


def resolve_block_node(row: object, lookup: dict) -> str:
    county, orig = str(row.COUNTY), str(row.NAME)
    sd, hd = int(row.ENACTED_UPPER), int(row.ENACTED_LOWER)
    key = (county, orig, sd, hd)
    if key in lookup:
        return lookup[key]
    # A cleaned block row can itself carry comma-joined precursor precincts.
    tokens = [piece.strip() for piece in orig.split(",")]
    if not all((county, piece, sd, hd) in lookup for piece in tokens):
        raise ValueError(f"Block precinct has unmapped components: {key}")
    matches = {lookup[(county, piece, sd, hd)] for piece in tokens}
    if len(matches) != 1:
        raise ValueError(f"Block precinct cannot be mapped uniquely: {key} -> {matches}")
    return next(iter(matches))


def spatial_membership(blocks: gpd.GeoDataFrame, selected: gpd.GeoDataFrame,
                       label: str) -> pd.Series:
    """At most one unit of each kind may contain a block representative point."""
    points = gpd.GeoDataFrame(geometry=blocks.geometry.representative_point(),
                              crs=blocks.crs, index=blocks.index)
    hits = gpd.sjoin(points, selected[["GEOID", "geometry"]],
                     how="left", predicate="within")
    if hits.index.duplicated().any():
        examples = hits.index[hits.index.duplicated()].tolist()[:5]
        raise ValueError(f"Blocks matched multiple {label} units: {examples}")
    return hits["GEOID"].reindex(blocks.index).fillna("").astype(str)


def build_atoms(blocks: gpd.GeoDataFrame, graph_nodes: list[dict],
                places: gpd.GeoDataFrame, townships: gpd.GeoDataFrame,
                classes: pd.DataFrame) -> tuple[list[dict], dict]:
    required = {"NAME", "COUNTY", "POP20", "ENACTED_UPPER", "ENACTED_LOWER", "geometry"}
    if not required.issubset(blocks.columns):
        raise ValueError(f"Block file missing columns {sorted(required - set(blocks.columns))}")
    if blocks.crs is None:
        raise ValueError("Blocks have no CRS")
    if not blocks.index.is_unique:
        raise ValueError("Block index must be unique")
    names, node_pop = graph_node_lookup(graph_nodes)
    # All geometry layers must match for representative-point membership.
    places = places.to_crs(blocks.crs)
    townships = townships.to_crs(blocks.crs)
    town_ids = spatial_membership(blocks, townships, "township")
    place_ids = spatial_membership(blocks, places, "incorporated place")
    place_status = classes.set_index("place_geoid")["status"].to_dict()
    aggregated = defaultdict(int)
    per_node = defaultdict(int)
    count = 0
    for i, row in enumerate(blocks.itertuples(index=False)):
        node = resolve_block_node(row, names)
        p = int(row.POP20)
        if p < 0:
            raise ValueError(f"Negative block population on {node}")
        per_node[node] += p
        t = town_ids.iat[i]
        m = place_ids.iat[i]
        classification = place_status[m] if m else "none"
        if p:
            aggregated[(node, t, m, classification)] += p
        count += 1
    wrong = [(node, actual, pop) for node, pop in node_pop.items()
             if (actual := per_node[node]) != pop]
    if wrong:
        raise ValueError(f"Population reconciliation failed for {len(wrong)} graph nodes: {wrong[:8]}")
    if sum(per_node.values()) != sum(node_pop.values()):
        raise AssertionError("Statewide graph and block populations disagree")
    atoms = [dict(zip(ATOM_FIELDS, (*key, pop)))
             for key, pop in sorted(aggregated.items())]
    stat = {
        "graph_nodes": len(node_pop),
        "census_block_rows": count,
        "aggregated_atoms": len(atoms),
        "statewide_pop20": sum(node_pop.values()),
        "township_pop20": sum(row["population"] for row in atoms if row["township_geoid"]),
        "nested_place_pop20": sum(row["population"] for row in atoms if row["place_class"] == "nested"),
        "cross_place_pop20": sum(row["population"] for row in atoms if row["place_class"] == "cross"),
        "n_townships": len(townships),
        "n_incorporated_places": len(places),
        "n_nested_places": int((classes.status == "nested").sum()),
        "n_cross_places": int((classes.status == "cross").sum()),
        "unassigned_township_pop20": sum(row["population"] for row in atoms if not row["township_geoid"]),
    }
    if stat["nested_place_pop20"] + stat["cross_place_pop20"] > stat["statewide_pop20"]:
        raise ValueError("Incorporated PLACE population exceeds statewide population")
    return atoms, stat


def parse_cli() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--blocks", default="data/derived/cblocksMergeMultiPcts_OH.gpkg")
    ap.add_argument("--graph", default="data/graph/OH_pct20.json")
    ap.add_argument("--places", default="data/raw/tl_2024_39_place")
    ap.add_argument("--cousubs", default="data/raw/tl_2024_39_cousub.zip")
    ap.add_argument("--atoms-out", required=True, help="TSV; one row per graph-node / unit intersection")
    ap.add_argument("--units-out", required=True, help="CSV audit of included geographic units and nesting")
    ap.add_argument("--summary-out", required=True, help="JSON population-reconciliation summary")
    ap.add_argument("--relative-sliver", type=float, default=1e-8)
    ap.add_argument("--absolute-sliver-m2", type=float, default=1.0)
    ap.add_argument("--overwrite", action="store_true")
    return ap.parse_args()


def main() -> None:
    args = parse_cli()
    out_paths = [Path(args.atoms_out), Path(args.units_out), Path(args.summary_out)]
    if len(set(out_paths)) != len(out_paths):
        raise ValueError("All output paths must be distinct")
    for path in out_paths:
        if path.exists() and not args.overwrite:
            raise FileExistsError(f"Refusing to overwrite {path}; pass --overwrite deliberately")
    places = selected_geographies(args.places, kind="place")
    townships = selected_geographies(args.cousubs, kind="township")
    places = places.to_crs("EPSG:26917")
    townships = townships.to_crs("EPSG:26917")
    classes = place_nesting(places, townships, relative_sliver=args.relative_sliver,
                            absolute_sliver_m2=args.absolute_sliver_m2)
    blocks = gpd.read_file(args.blocks)
    graph = json.loads(Path(args.graph).read_text())
    atoms, stats = build_atoms(blocks, graph["nodes"], places, townships, classes)
    units = []
    for t in townships.itertuples(index=False):
        units.append({"geography": "township", "geoid": t.GEOID, "name": t.NAME,
                      "classfp": t.CLASSFP, "lsad": t.LSAD, "status": "township",
                      "parent_township_geoid": "", "nearest_township_geoid": "",
                      "outside_nearest_township_area_ratio": ""})
    for c in classes.itertuples(index=False):
        units.append({"geography": "place", "geoid": c.place_geoid,
                      "name": c.place_name, "classfp": c.place_classfp, "lsad": c.place_lsad,
                      "status": c.status, "parent_township_geoid": c.parent_township_geoid,
                      "nearest_township_geoid": c.nearest_township_geoid,
                      "outside_nearest_township_area_ratio": c.outside_nearest_township_area_ratio})
    for path in out_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    # Commit data files only after all validations pass. No partial successful CSVs.
    temp_paths = [Path(str(path) + ".tmp") for path in out_paths]
    try:
        with temp_paths[0].open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=ATOM_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(atoms)
        with temp_paths[1].open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=("geography", "geoid", "name", "classfp", "lsad", "status",
                            "parent_township_geoid", "nearest_township_geoid", "outside_nearest_township_area_ratio"))
            writer.writeheader()
            writer.writerows(sorted(units, key=lambda u: (u["geography"], u["geoid"])))
        temp_paths[2].write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
        for temp, target in zip(temp_paths, out_paths):
            temp.replace(target)
    finally:
        for path in temp_paths:
            path.unlink(missing_ok=True)
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
