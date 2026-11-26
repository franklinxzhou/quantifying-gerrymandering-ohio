# Record of Commands

The archived launchers specify all four treatments and their output paths. The recovered Julia source agrees with the dated implementation bundles described below. The archived manifest identifies the exact dependency source tree available in the recovered checkout. 

## Original commands and chain locations

The commands below retain the historical argument lists. They ran from the original `oh_gerrymandering` directory on Duke Math Cluster. In these excerpts, `OUT` is the chain directory in the table, relative to the authors' scratch directory. The statewide commands used `BASE` for the `oh-mct-pilot-20261118` scratch directory and `GEOM="$BASE/mct_edge_geometry_v3"`. The corresponding portable commands are in `sampling/`.

| Ensemble | Chain directory under the scratch root | Filenames |
| --- | --- | --- |
| Minimal | `oh-congressional-columbus-county-energy-20261118/production/lambda-5.0/chains/freq-100` | `columbus_county_lambda5.0_run{8601,8602,8603,8604}_cntyw50_placew5_B7.jsonl` |
| Match 2025 | `oh-mct-production-mc10-t1-l10-20261120/chains` | `mc10_t1_l10-run{1,2,3,4}.jsonl` |
| Match 2021/22 | `oh-mct-production-mc10-t5-l10-20261120/chains` | `mc10_t5_l10-run{1,2,3,4}.jsonl` |
| High Preservation | `oh-mct-production-mc20-t5-l20-lt5-20261121/chains` | `mc20_t5_l20_lt5-run{1,2,3,4}.jsonl` |

### Minimal

The launcher set `JULIA_NUM_THREADS=1` and `JULIA_PKG_PRECOMPILE_AUTO=0`. It used Sun Grid Engine task IDs 8601–8604 and verified the 5,000-map quota after each run.

```bash
ID="$SGE_TASK_ID"
FILE="columbus_county_lambda5.0_run${ID}_cntyw50_placew5_B7.jsonl"
julia --project=. scripts/grid/columbus_baseline_county_energy/lct_g0_columbus.jl \
  --columbus_only --chamber congressional --preserve_five_cities --place_county_units \
  --enforce_b7 --place_geoids_file scripts/grid/columbus_baseline_county_energy/columbus_geoid.txt \
  --place_w 5.0 --cnty_w 50.0 --place_ousted_weight 5.0 --steps 20000000 \
  --output_freq 100 --target_maps 5000 --max_districts_per_county 2 --popdev 0.02 \
  --gamma 0.0 --iso_weight 0.0 --const_weight 0.0 \
  --outpath "$OUT" --outfile "$FILE" --thread "$ID"
```

### Match 2025

The production analysis export points to these filenames. 

```bash
THREAD=$((120 + CHAIN))
FILE="mc10_t1_l10-run${CHAIN}.jsonl"
julia --project=. lct_g0_mct.jl \
  --chamber congressional --target_maps 5000 --output_freq 100 \
  --thread "$THREAD" --cnty_w 50 --place_w 1 --popdev 0.02 \
  --mct_atoms_file "$BASE/mct_pop_nested_atoms.tsv" --mct_edge_model full \
  --mct_cross_place_w 10 --mct_township_w 1 \
  --mct_edge_cross_ids_file "$GEOM/geometric_cross_places.txt" \
  --mct_edge_incidence_file "$GEOM/node_geom_incidence.tsv" \
  --mct_cross_place_ousted_lambda 10 --mct_township_ousted_lambda 0 \
  --max_districts_per_county 2 --outpath "$OUT" --outfile "$FILE"
```

### Match 2021/22

```bash
THREAD=$((110 + CHAIN))
FILE="mc10_t5_l10-run${CHAIN}.jsonl"
julia --project=. lct_g0_mct.jl \
  --chamber congressional --target_maps 5000 --output_freq 100 \
  --thread "$THREAD" --cnty_w 50 --place_w 1 --popdev 0.02 \
  --mct_atoms_file "$BASE/mct_pop_nested_atoms.tsv" --mct_edge_model full \
  --mct_cross_place_w 10 --mct_township_w 5 \
  --mct_edge_cross_ids_file "$GEOM/geometric_cross_places.txt" \
  --mct_edge_incidence_file "$GEOM/node_geom_incidence.tsv" \
  --mct_cross_place_ousted_lambda 10 --mct_township_ousted_lambda 0 \
  --max_districts_per_county 2 --outpath "$OUT" --outfile "$FILE"
```

### High Preservation

```bash
THREAD=$((130 + CHAIN))
FILE="mc20_t5_l20_lt5-run${CHAIN}.jsonl"
julia --project=. lct_g0_mct.jl \
  --chamber congressional --target_maps 5000 --output_freq 100 \
  --thread "$THREAD" --cnty_w 50 --place_w 1 --popdev 0.02 \
  --mct_atoms_file "$BASE/mct_pop_nested_atoms.tsv" --mct_edge_model full \
  --mct_cross_place_w 20 --mct_township_w 5 \
  --mct_edge_cross_ids_file "$GEOM/geometric_cross_places.txt" \
  --mct_edge_incidence_file "$GEOM/node_geom_incidence.tsv" \
  --mct_cross_place_ousted_lambda 20 --mct_township_ousted_lambda 5 \
  --max_districts_per_county 2 --outpath "$OUT" --outfile "$FILE"
```

### Defaults and seeds

Both recovered drivers use `PCG.PCGStateOneseq(UInt64, 1029472 + 512339 * thread)`. The seed is used for the random initial partition, the weighted spanning trees that initialize its districts, and the subsequent chain. CycleWalk's constructor first draws a population-balanced `MultiLevelPartition` using MetropolizedForestRecom, then converts the assignment into a `LinkCutPartition` with trees drawn using Wilson's algorithm. Neither recovered launcher supplies an initial assignment file. [Pinned CycleWalk initialization](https://github.com/jonmjonm/CycleWalk.jl/blob/f24f129f17d78707fd05931e83a298c54cd598da/src/partition/link_cut_partition.jl)

| Ensemble | Values of `--thread` | RNG seeds in chain order |
| --- | --- | --- |
| Minimal | 8601–8604 | 4407657211, 4408169550, 4408681889, 4409194228 |
| Match 2025 | 121–124 | 63022491, 63534830, 64047169, 64559508 |
| Match 2021/22 | 111–114 | 57899101, 58411440, 58923779, 59436118 |
| High Preservation | 131–134 | 68145881, 68658220, 69170559, 69682898 |

The statewide launch commands rely on defaults of 40,000,000 nominal steps, `gamma=0`, `iso_weight=0`, and `const_weight=5`. The constitutional energy callback for `const_weight` is specific to the state upper chamber. Minimal explicitly sets the nominal step argument to 20,000,000 and all three of those coefficients to zero. Its county-subdivision energy defaults to zero. All four drivers stop at the retained-map quota. 

