# Ohio congressional redistricting ensembles

This repository accompanies *Quantifying Gerrymandering in Ohio: A Cycle Walk-Based Ensemble Analysis* by G.J. Herschlag, J.C. Mattingly, and F.X. Zhou. The [record of commands](PROVENANCE.md) gives the original launch commands.

## Installation

Use Julia 1.12.1, the version recorded in the recovered manifest, and Bash. From the repository root, install the pinned environment and check the sampler dependencies.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using CycleWalk, RandomNumbers'
```

The manifest pins CycleWalk v0.2.0. Its source tree, `6d9d49eeb8680513178be645231fcc497c02b6b2`, matches upstream commit `f24f129f17d78707fd05931e83a298c54cd598da`. The recovered Grid environment and the production jobs are documented separately in the [production record](PROVENANCE.md). [CycleWalk release commit](https://github.com/jonmjonm/CycleWalk.jl/commit/f24f129f17d78707fd05931e83a298c54cd598da)

## Geographic inputs

- `data/graph/OH_pct20.json` is the recovered Ohio adjacency graph with 2020 Census population (`POP20`), 7,928 precinct units after articulation merges, and 11,799,448 residents. The statewide ensembles use this graph directly. Its `NAME` keys identify precincts in all four Atlas ensembles. [2020 Census redistricting data](https://www.census.gov/programs-surveys/decennial-census/about/rdo/summary-files.html)
- `data/graph/OH_pct20_preserved_five_cities_columbus_only.json` is the Minimal sampling graph. It has 6,975 nodes, including one contracted node each for Akron, Cincinnati, Cleveland, Dayton, and Toledo. `OH_pct20_preserved_five_cities.json` supplies the original contraction and its `MEMBER_NAMES` mapping for rebuilding the Minimal graph. `columbus_geoid.txt` selects Columbus (PLACE GEOID `3918000`).
- `data/mct/` contains the population-effective municipality and township atoms, their classification and population summaries, the geometrically non-nested PLACE IDs, and the node incidence table. The statewide samplers use the prepared atoms and the two geometric files.

### Rebuild the geographic inputs

The prepared inputs above suffice for sampling. To regenerate the Minimal graph, run the archived PLACE-by-county enrichment and Columbus filtering on the packaged graphs. The enrichment verifies the existing contraction by mapping every original precinct to a retained precinct or protected-city node and checking population and edge topology.

```bash
mkdir -p data/rebuilt/graph
python3 preprocessing/build_oh_place_county_graph.py \
  --original-output data/rebuilt/graph/OH_pct20_place_county.json \
  --contracted-output data/rebuilt/graph/OH_pct20_preserved_five_cities_place_county.json \
  --report data/rebuilt/graph/place-county-audit.json
python3 preprocessing/build_columbus_graph.py \
  --precinct data/rebuilt/graph/OH_pct20_place_county.json \
  --five-city data/rebuilt/graph/OH_pct20_preserved_five_cities_place_county.json \
  --output-dir data/rebuilt/graph
```

To rebuild the statewide MCT inputs, use Python 3.9 or newer, the Ohio 2024 TIGER/Line PLACE and COUSUB layers, and the historical block-to-precinct GeoPackage `cblocksMergeMultiPcts_OH.gpkg`. Obtain the GeoPackage from the authors' Grid archive and place it in `data/derived/`. It contains merged Census block records with geometry, original precinct name, county, POP20, and enacted state House and Senate district labels. The atoms builder uses those fields to resolve the final precinct name and checks population against every graph node. [2024 TIGER/Line files](https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.2024.html)

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r preprocessing/requirements.txt
mkdir -p data/raw data/derived data/rebuilt/mct

curl -fLo data/raw/tl_2024_39_place.zip \
  https://www2.census.gov/geo/tiger/TIGER2024/PLACE/tl_2024_39_place.zip
curl -fLo data/raw/tl_2024_39_cousub.zip \
  https://www2.census.gov/geo/tiger/TIGER2024/COUSUB/tl_2024_39_cousub.zip
unzip -n data/raw/tl_2024_39_place.zip -d data/raw/tl_2024_39_place
unzip -n data/raw/tl_2024_39_cousub.zip -d data/raw/tl_2024_39_cousub

graph=data/graph/OH_pct20.json
blocks=data/derived/cblocksMergeMultiPcts_OH.gpkg
places=data/raw/tl_2024_39_place/tl_2024_39_place.shp
townships=data/raw/tl_2024_39_cousub/tl_2024_39_cousub.shp
out=data/rebuilt/mct

python preprocessing/build_mct_atoms.py \
  --graph "$graph" --blocks "$blocks" --places "$places" --cousubs "$townships" \
  --atoms-out "$out/mct_atoms.tsv" --units-out "$out/mct_units.csv" \
  --summary-out "$out/mct_summary.json"
python preprocessing/mct_population_nesting.py \
  --atoms "$out/mct_atoms.tsv" --units "$out/mct_units.csv" \
  --summary "$out/mct_summary.json" --output-prefix "$out/mct_pop_nested"
python preprocessing/build_edge_geometry_incidence.py \
  --graph "$graph" --blocks "$blocks" --places "$places" --townships "$townships" \
  --units "$out/mct_pop_nested_units.csv" \
  --cross-ids-out "$out/geometric_cross_places.txt" \
  --incidence-out "$out/node_geom_incidence.tsv"
```

To sample with rebuilt inputs, compare them with the prepared files and use their paths in the Julia invocations in `sampling/`. The historical upstream block and precinct construction and five-city contraction are described in the [record of commands](PROVENANCE.md).