#!/usr/bin/env bash
set -euo pipefail

# SGE spools the script. Submit array jobs from the repository root.
if [[ -n ${JOB_ID:-} ]]; then
    cd "${SGE_O_WORKDIR:?Submit from the repository root with qsub -cwd}"
else
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
fi

run=${1:-${SGE_TASK_ID:-}}
case "$run" in
    8601|8602|8603|8604) ;;
    *) echo "Usage: $0 RUN_ID (8601 through 8604)" >&2; exit 2 ;;
esac

out="${OUTPUT_DIR:-outputs}/alpha5"
file="columbus_county_lambda5.0_run${run}_cntyw50_placew5_B7.jsonl"
exec "${JULIA:-julia}" --project=. --threads=1 sampling/lct_g0_columbus.jl \
    --chamber congressional --columbus_only --preserve_five_cities \
    --place_county_units --enforce_b7 --place_geoids_file data/graph/columbus_geoid.txt \
    --place_w 5.0 --cnty_w 50.0 --place_ousted_weight 5.0 --cousub_ousted_weight 0.0 \
    --popdev 0.02 --max_districts_per_county 2 --gamma 0.0 --iso_weight 0.0 --const_weight 0.0 \
    --steps 20000000 --output_freq 100 --target_maps 5000 \
    --base_seed 1029472 --mul_seed 512339 --thread "$run" \
    --outpath "$out" --outfile "$file"
