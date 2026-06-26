#!/bin/bash
#SBATCH --job-name=3way_tmix_agg
#SBATCH --time=01:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A.err

## Walk data/<MODE>/tractormix_out/alt/* and join with SAIGE outputs.
##   MODE=common  sbatch 18_aggregate_tractormix.sh

set -euo pipefail
MODE=${MODE:?MODE must be exported}
module load singularity

BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/scripts/R/aggregate_tractormix.R "${MODE}"
$SING $RTOOLS Rscript ${BASE}/scripts/R/plot_tractormix_bench.R  "${MODE}"
