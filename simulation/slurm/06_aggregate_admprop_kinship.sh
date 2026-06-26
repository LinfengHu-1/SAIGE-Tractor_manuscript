#!/bin/bash
#SBATCH --job-name=3way_agg_adm
#SBATCH --time=01:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A.err

## Aggregate per-block Admprop files into Admprop.tsv AND write the true
## pedigree kinship matrix to Kinship.tsv (used by phenotype generators).
##
##   MODE=common  sbatch 06_aggregate_admprop_kinship.sh
##   MODE=lowfreq sbatch 06_aggregate_admprop_kinship.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/scripts/R/aggregate_Admprop.R "${MODE}"
