#!/bin/bash -e
#SBATCH --job-name=Yale25_trimmomatic_alignment_array
#SBATCH --time=7-00:00:00
#SBATCH --array=1-18
#SBATCH --output=/work/hs325/bass25/2024/misc/Yale25_trimmomatic_hisat2_alignment_%a.out
#SBATCH --error=/work/hs325/bass25/2024/misc/Yale25_trimmomatic_hisat2_alignment_%a.err
#SBATCH --partition=schultzlab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

# source conda
source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh

## Load module - HISAT2 already loaded on DCC ## 
module load HISAT2

## Set Paths for Alignment ## 
RAW_DIR=/work/hs325/bass25/2024/raw
TRIMMED_DIR=/work/hs325/bass25/2024/trimmed
INDEX_DIR=/work/hs325/bass25/align/index
ALIGNED_DIR=/work/hs325/bass25/2024/aligned
mkdir -p ${ALIGNED_DIR} 

## Set up direction/path to each sample ##
# Make list of trimmed sample names
SAMPLES=($(ls ${RAW_DIR}/*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' | xargs -n 1 basename))

# Index an individual sample from the list for this array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}
# List the sample to see if naming the correct thing 
echo "Sample being processed here: " ${SAMPLE}

# Define R1 and R2 for the sample 
R1=${TRIMMED_DIR}/${SAMPLE}_R1.paired.fastq.gz
R2=${TRIMMED_DIR}/${SAMPLE}_R2.paired.fastq.gz
# List R1 and R2 to see if naming the correct thing 
echo "Path to R1 is " ${R1} "and path to R2 is " ${R2}

## Run Alignment ##
echo "Aligning sample:" ${SAMPLE}
hisat2 \
    -p ${SLURM_CPUS_PER_TASK} \
    -x ${INDEX_DIR}/c.virginica_yale25_HFM_index \
    -1 ${R1} \
    -2 ${R2} \
    -S ${ALIGNED_DIR}/${SAMPLE}.sam \
    --summary-file ${ALIGNED_DIR}/${SAMPLE}_hisat2_summary.txt

echo "Alignment of " ${SAMPLE} "complete!"
