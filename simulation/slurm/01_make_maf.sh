#!/bin/bash
#SBATCH --job-name=3way_maf
#SBATCH --time=01:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A.err

## Generate the AFR/EUR/NAT ancestral and population-specific allele frequencies.
## Submit one job per MAF mode:
##   MODE=common  sbatch 01_make_maf.sh
##   MODE=lowfreq sbatch 01_make_maf.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity
BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/scripts/R/Pedigree3way_MAF.R "${MODE}"
