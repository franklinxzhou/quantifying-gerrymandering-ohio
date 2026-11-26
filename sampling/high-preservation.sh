#!/usr/bin/env bash
set -euo pipefail

# SGE spools the script. Submit array jobs from the repository root.
if [[ -n ${JOB_ID:-} ]]; then
    cd "${SGE_O_WORKDIR:?Submit from the repository root with qsub -cwd}"
else
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
fi

chain=${1:-${SGE_TASK_ID:-}}
case "$chain" in
    1|2|3|4) ;;
    *) echo "Usage: $0 CHAIN (1 through 4)" >&2; exit 2 ;;
esac

thread=$((130 + chain))
out="${OUTPUT_DIR:-outputs}/mc20_t5_l20_lt5"
file="mc20_t5_l20_lt5-run${chain}.jsonl"
exec "${JULIA:-julia}" --project=. --threads=1 sampling/lct_g0_mct.jl \
    --chamber congressional --cnty_w 50 --place_w 1 --popdev 0.02 --max_districts_per_county 2 \
    --mct_atoms_file data/mct/mct_pop_nested_atoms.tsv --mct_edge_model full \
    --mct_edge_cross_ids_file data/mct/geometric_cross_places.txt \
    --mct_edge_incidence_file data/mct/node_geom_incidence.tsv \
    --mct_cross_place_w 20 --mct_township_w 5 \
    --mct_cross_place_ousted_lambda 20 --mct_township_ousted_lambda 5 \
    --gamma 0.0 --iso_weight 0.0 --const_weight 5.0 \
    --steps 40000000 --output_freq 100 --target_maps 5000 \
    --base_seed 1029472 --mul_seed 512339 --thread "$thread" \
    --outpath "$out" --outfile "$file"
