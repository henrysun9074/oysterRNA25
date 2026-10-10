#!/bin/bash -e
#SBATCH --job-name=qc
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/hs325/bass25/misc/multiqc.out
#SBATCH --error=/work/hs325/bass25/misc/multiqc.err
#SBATCH --partition=scavenger
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate RNA-seq

# Load module
module load FastQC

## Set paths ##
RAW_DIR=/work/hs325/bass25/clean
FASTQC_OUT=/work/hs325/bass25/clean
mkdir -p $FASTQC_OUT

## Set up direction/path to each sample ##
# Make list of sample names (without _R1/_R2 suffix)
SAMPLES=($(ls ${RAW_DIR}/*_R1.trimmed.fastq.gz | sed 's/_R1.trimmedfastq.gz//' | xargs -n 1 basename))

# Index an individual sample from the list for this array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Define R1 and R2 for the sample 
R1=${RAW_DIR}/${SAMPLE}_R1.trimmed.fastq.gz
R2=${RAW_DIR}/${SAMPLE}_R2.trimmed.fastq.gz

## Run FastQC ##
echo "Running FastQC for sample: ${SAMPLE}"
fastqc -o ${FASTQC_OUT} -t ${SLURM_CPUS_PER_TASK} ${R1} ${R2}

echo "FastQC complete for sample: ${SAMPLE}"

###### run multiqc
MULTIQC_OUT=/work/hs325/bass25/clean/multiqc_trimmed
mkdir -p $MULTIQC_OUT

multiqc $FASTQC_OUT -o $MULTIQC_OUT --force
echo "MultiQC complete!"

conda deactivate