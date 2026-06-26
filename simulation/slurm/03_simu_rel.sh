#!/bin/bash
#SBATCH --job-name=3way_rel
#SBATCH --time=12:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.err
#SBATCH --array=1-50

## 50 array tasks x (10 families x 10 members) = 5000 related admixed samples.
##   MODE=common  sbatch 03_simu_rel.sh
##   MODE=lowfreq sbatch 03_simu_rel.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/scripts/R/Pedigree3way_Rel.R "${SLURM_ARRAY_TASK_ID}" "${MODE}"
