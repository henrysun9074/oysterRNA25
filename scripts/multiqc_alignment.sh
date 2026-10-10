#!/bin/bash
#SBATCH --job-name=postalign_multiqc
#SBATCH --time=02:00:00
#SBATCH --output=/work/hs325/bass25/misc/postalign_multiqc.out
#SBATCH --error=/work/hs325/bass25/misc/postalign_multiqc.err
#SBATCH --partition=scavenger
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

set -eo pipefail

source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate RNA-seq

SAM_DIR="${1:-/work/hs325/bass25/align}"
MULTIQC_OUT="$SAM_DIR/multiqc"

shopt -s nullglob
summaries=("$SAM_DIR"/*_hisat2_summary.txt)
if (( ${#summaries[@]} == 0 )); then
    echo "No HISAT2 summary files found in $SAM_DIR" >&2
    exit 1
fi

echo "Running MultiQC on ${#summaries[@]} HISAT2 summary files"
multiqc "${summaries[@]}" -o "$MULTIQC_OUT"
echo "MultiQC complete: $MULTIQC_OUT"
