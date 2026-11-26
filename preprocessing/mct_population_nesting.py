#!/usr/bin/env python3
"""Derive population-effective nesting from an existing validated Ohio MCT crosswalk.

Retains polygon-based classification as `geometric_status`; produces an
alternative atoms file consumable by the existing mct_energy.jl scorer.
Original files are never modified. No geospatial dependencies are required.
"""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ATOM_FIELDS = ["node_name", "township_geoid", "place_geoid", "place_class", "population"]
UNIT_FIELDS = ["geography", "geoid", "name", "classfp", "lsad", "status", "parent_township_geoid", "nearest_township_geoid", "outside_nearest_township_area_ratio"]


def read_csv(path, *, delimiter=","):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        expected = ATOM_FIELDS if delimiter == "\t" else UNIT_FIELDS
        if reader.fieldnames != expected:
            raise ValueError(f"Unexpected columns in {path}: {reader.fieldnames}")
        return list(reader)


def convert(atoms, units, original_summary):
    place_units = {r["geoid"]: r for r in units if r["geography"] == "place"}
    township_units = {r["geoid"] for r in units if r["geography"] == "township"}
    if len(place_units) != sum(r["geography"] == "place" for r in units):
        raise ValueError("Duplicate PLACE GEOID")
    if len(township_units) != sum(r["geography"] == "township" for r in units):
        raise ValueError("Duplicate township GEOID")
    totals, outside = Counter(), Counter()
    node_totals = Counter()
    atomkeys = set()
    for row in atoms:
        node, township, place, label = (row[x] for x in ATOM_FIELDS[:4])
        population = int(row["population"])
        if population <= 0:
            raise ValueError(f"Nonpositive atom population: {row}")
        key = (node, township, place)
        if key in atomkeys:
            raise ValueError(f"Duplicate population atom {key}")
        atomkeys.add(key)
        node_totals[node] += population
        if township and township not in township_units:
            raise ValueError(f"Unknown township {township}")
        if not place:
            if label != "none":
                raise ValueError(f"Missing place with class {label}")
            continue
        unit = place_units.get(place)
        if unit is None or label != unit["status"]:
            raise ValueError(f"Inconsistent atom and unit classification: {place} {label}")
        totals[place] += population
        parent = (unit["parent_township_geoid"] if label == "nested"
                  else unit["nearest_township_geoid"])
        if label == "nested" and parent not in township_units:
            raise ValueError(f"Nested place {place} lacks qualifying parent township")
        if not parent or township != parent:
            outside[place] += population
    if sum(node_totals.values()) != int(original_summary["statewide_pop20"]):
        raise ValueError("Statewide population differs from validated summary")
    if len(node_totals) > int(original_summary["graph_nodes"]):
        raise ValueError("More populated nodes than original graph nodes")
    if len(atoms) != int(original_summary["aggregated_atoms"]):
        raise ValueError("Original atom count differs from validated summary")

    promoted = set()
    for geoid, unit in place_units.items():
        status = unit["status"]
        if status not in ("nested", "cross"):
            raise ValueError(f"Unexpected place status: {status}")
        if status == "nested" and outside[geoid]:
            raise ValueError(f"Geometrically nested {geoid} has population outside its parent")
        if status == "cross" and totals[geoid] > 0 and outside[geoid] == 0:
            candidate = unit["nearest_township_geoid"]
            if candidate not in township_units:
                raise ValueError(f"Proposed promoted place {geoid} lacks a valid township")
            promoted.add(geoid)

    transformed_atoms = []
    for row in atoms:
        new = row.copy()
        if row["place_geoid"] in promoted:
            new["place_class"] = "nested"
        transformed_atoms.append(new)

    transformed_units = []
    for unit in units:
        new = unit.copy()
        if unit["geography"] == "place":
            new["geometric_status"] = unit["status"]
            if unit["status"] == "nested":
                new["nesting_basis"] = "geometry"
            elif unit["geoid"] in promoted:
                new["status"] = "nested"
                new["parent_township_geoid"] = new["nearest_township_geoid"]
                new["nesting_basis"] = "population"
            else:
                new["nesting_basis"] = "non_nested"
        else:
            new["geometric_status"] = ""
            new["nesting_basis"] = ""
        transformed_units.append(new)

    old_nested = sum(u["status"] == "nested" for u in place_units.values())
    new_nested = old_nested + len(promoted)
    summary = original_summary.copy()
    summary.update({
        "nesting_mode": "population_effective",
        "n_geometry_nested_places": old_nested,
        "n_population_promoted_places": len(promoted),
        "n_nested_places": new_nested,
        "n_cross_places": len(place_units) - new_nested,
        "n_nested_places_with_population": sum(totals[p] > 0 for p, u in place_units.items() if u["status"] == "nested") + len(promoted),
        "nested_place_pop20": int(sum(totals[p] for p,u in place_units.items() if u["status"] == "nested" or p in promoted)),
        "cross_place_pop20": int(sum(totals[p] for p,u in place_units.items() if u["status"] == "cross" and p not in promoted)),
        "promoted_place_pop20": int(sum(totals[p] for p in promoted)),
    })
    if summary["nested_place_pop20"] + summary["cross_place_pop20"] != (int(original_summary["nested_place_pop20"]) + int(original_summary["cross_place_pop20"])):
        raise ValueError("PLACE population not conserved by promotion")
    return transformed_atoms, transformed_units, summary, sorted(promoted)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--atoms", required=True)
    ap.add_argument("--units", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--output-prefix", required=True)
    args = ap.parse_args()
    prefix = Path(args.output_prefix)
    output_paths = [Path(str(prefix) + suffix) for suffix in ("_atoms.tsv", "_units.csv", "_summary.json", "_promoted.txt")]
    for path in output_paths:
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite existing {path}")
    atoms = read_csv(args.atoms, delimiter="\t")
    units = read_csv(args.units)
    summary = json.loads(Path(args.summary).read_text())
    revised_atoms, revised_units, revised_summary, promoted = convert(atoms, units, summary)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    try:
        with output_paths[0].open("x", newline="") as f:
            w = csv.DictWriter(f, fieldnames=ATOM_FIELDS, delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(revised_atoms)
        with output_paths[1].open("x", newline="") as f:
            w = csv.DictWriter(f, fieldnames=UNIT_FIELDS + ["geometric_status", "nesting_basis"], lineterminator="\n")
            w.writeheader(); w.writerows(revised_units)
        with output_paths[2].open("x") as f:
            json.dump(revised_summary, f, indent=2, sort_keys=True)
            f.write("\n")
        output_paths[3].write_text("\n".join(promoted) + "\n")
    except Exception:
        for path in output_paths:
            path.unlink(missing_ok=True)
        raise
    print(json.dumps(revised_summary, indent=2, sort_keys=True))
    print(f"Promotion report: {output_paths[3]}")

if __name__ == "__main__":
    main()
