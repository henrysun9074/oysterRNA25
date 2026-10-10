#!/bin/bash 
#SBATCH --job-name=rasflow
#SBATCH --time=7-00:00:00
#SBATCH --chdir=/work/hs325/bass25/rasflow
#SBATCH --output=/work/hs325/bass25/misc/rasflow.out
#SBATCH --error=/work/hs325/bass25/misc/rasflow.err
#SBATCH --partition=schultzlab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

# set -euo pipefail

source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate rasflow

run_stage() {
    local stage="$1"
    snakemake \
        --snakefile "workflow/${stage}.rules" \
        --cores 1 \
        --printshellcmds \
        2>&1 | tee "logs/${stage}-${SLURM_JOB_ID}.log"
}

run_stage quality_control    || exit "$?"
run_stage trim               || exit "$?"
run_stage align_count_genome  || exit "$?"
run_stage dea_genome          || exit "$?"
run_stage visualize          || exit "$?"