#!/bin/bash
#SBATCH --job-name=postalign
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/hs325/bass25/misc/postalign_%a.out
#SBATCH --error=/work/hs325/bass25/misc/postalign_%a.err
#SBATCH --partition=scavenger
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

set -eo pipefail

ALIGNED_DIR="/work/hs325/bass25/align"
BAM_DIR="$ALIGNED_DIR/bam"
LOG_DIR="/work/hs325/bass25/misc"
SORT_MEMORY_PER_THREAD="2G"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    shopt -s nullglob
    sam_files=("$ALIGNED_DIR"/*.sam)
    n_samples="${#sam_files[@]}"
    if (( n_samples == 0 )); then
        echo "No SAM files found in $ALIGNED_DIR" >&2
        exit 1
    fi

    mkdir -p "$BAM_DIR" "$LOG_DIR"
    manifest="$(mktemp "$LOG_DIR/postalign_samples.XXXXXX")"
    printf '%s\n' "${sam_files[@]}" > "$manifest"
    script_path="$(readlink -f "${BASH_SOURCE[0]}")"

    submission="$(sbatch --parsable \
        --array="1-${n_samples}" \
        --export="ALL,POSTALIGN_MANIFEST=$manifest" \
        "$script_path")"
    echo "Submitted postalign array ${submission%%;*} with $n_samples samples"
    echo "SAM-file list: $manifest"
    exit 0
fi

manifest="${POSTALIGN_MANIFEST:?Run this workflow with bash postalign_array.sh}"
task_id="${SLURM_ARRAY_TASK_ID:?This script requires a SLURM array}"
threads="${SLURM_CPUS_PER_TASK:-8}"

mapfile -t sam_files < "$manifest"
if (( task_id < 1 || task_id > ${#sam_files[@]} )); then
    echo "Invalid array index: $task_id" >&2
    exit 1
fi
sam="${sam_files[$((task_id - 1))]}"
SAMPLE="$(basename "$sam" .sam)"
if [[ ! -s "$sam" ]]; then
    echo "Missing or empty SAM: $sam" >&2
    exit 1
fi

module load samtools
mkdir -p "$BAM_DIR"
unsorted_bam="$BAM_DIR/${SAMPLE}.bam"
sorted_bam="$BAM_DIR/${SAMPLE}_sorted.bam"

echo "[$SAMPLE] Converting SAM to BAM"
samtools view -b "$sam" -o "$unsorted_bam"

echo "[$SAMPLE] Sorting BAM"
samtools sort -@ "$threads" -m "$SORT_MEMORY_PER_THREAD" \
    "$unsorted_bam" -o "$sorted_bam"

samtools quickcheck -v "$sorted_bam"
samtools index "$sorted_bam"
echo "[$SAMPLE] Finished BAM conversion, sorting and indexing"
