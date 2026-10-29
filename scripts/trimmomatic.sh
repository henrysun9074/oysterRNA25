#!/bin/bash -e
#SBATCH --job-name=trimmomatic
#SBATCH --time=7-00:00:00
#SBATCH --array=1-27
#SBATCH --output=/work/hs325/bass25/misc/trimmomatic_%a.out
#SBATCH --error=/work/hs325/bass25/misc/trimmomatic_%a.err
#SBATCH --partition=scavenger
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## Load fastp (poly-G trimming only) and Trimmomatic ##
source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate RNA-seq
module load Trimmomatic/0.38

## Set paths ##
RAW_DIR=/work/hs325/bass25/raw
TRIMMED_DIR=/work/hs325/bass25/clean/trimmomatic
ADAPTERS=/opt/apps/rhel7/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa
mkdir -p ${TRIMMED_DIR}

## Set up direction/path to each sample ##
# Make list of raw sample names
SAMPLES=($(ls ${RAW_DIR}/*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' | xargs -n 1 basename))

# Index an individual sample from the list for this array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}
echo "Processing sample: " ${SAMPLE}

# Define R1 and R2 for the sample
R1=${RAW_DIR}/${SAMPLE}_R1_001.fastq.gz
R2=${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz
echo "Path to R1 is " ${R1} "and path to R2 is " ${R2}

## Remove poly-G tails with fastp 
POLYG_R1=${TRIMMED_DIR}/${SAMPLE}_R1.polyG.fastq.gz
POLYG_R2=${TRIMMED_DIR}/${SAMPLE}_R2.polyG.fastq.gz
fastp -i ${R1} -I ${R2} \
      -o ${POLYG_R1} \
      -O ${POLYG_R2} \
      --trim_poly_g \
      --poly_g_min_len 10 \
      --disable_adapter_trimming \
      --disable_quality_filtering \
      --disable_length_filtering \
      -h ${TRIMMED_DIR}/${SAMPLE}.polyG.fastp.html \
      -j ${TRIMMED_DIR}/${SAMPLE}.polyG.fastp.json \
      --thread ${SLURM_CPUS_PER_TASK}

## Run Trimmomatic for paired-end data ##
java -Xmx12g -jar ${TRIMMOMATIC} PE \
    -threads ${SLURM_CPUS_PER_TASK} \
    -phred33 \
    -summary ${TRIMMED_DIR}/${SAMPLE}.trimmomatic_summary.txt \
    ${POLYG_R1} ${POLYG_R2} \
    ${TRIMMED_DIR}/${SAMPLE}_R1.paired.fastq.gz \
    ${TRIMMED_DIR}/${SAMPLE}_R1.unpaired.fastq.gz \
    ${TRIMMED_DIR}/${SAMPLE}_R2.paired.fastq.gz \
    ${TRIMMED_DIR}/${SAMPLE}_R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10:2:True \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:30 \
    2> ${TRIMMED_DIR}/${SAMPLE}.trimmomatic.log

# Remove the intermediate poly-G-trimmed files
rm ${POLYG_R1} ${POLYG_R2}

echo "Trimming of " ${SAMPLE} "complete!"
