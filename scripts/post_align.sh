#!/bin/bash -e
#SBATCH --job-name=postalign
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/hs325/bass25/misc/postalign.out
#SBATCH --error=/work/hs325/bass25/misc/postalign.err
#SBATCH --partition=scavenger
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

######## sam to bam ##########
module load samtools
echo samtools

## Set Paths ##
ALIGNED_DIR=/work/hs325/bass25/align
mkdir -p ${ALIGNED_DIR} 
SAM_DIR=${ALIGNED_DIR}
BAM_DIR=${SAM_DIR}/bam/
mkdir -p ${BAM_DIR}

## Set up direction/path to each sample ##
# Make list of trimmed sample names (without _R1/_R2 suffix)
SAMPLES=($(ls ${SAM_DIR}/*.sam | sed 's/.sam//' | xargs -n 1 basename))

# Index an individual sample from the list for this array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

## Convert from SAM to BAM file format ##
echo "Converting sample " ${SAMPLE} " from SAM to BAM..."

# First convert SAM --> BAM
samtools view -b ${SAM_DIR}/${SAMPLE}.sam -o ${BAM_DIR}/${SAMPLE}.bam
# Next sort the alignment (BAM) file
samtools sort -@ ${SLURM_CPUS_PER_TASK} ${BAM_DIR}/${SAMPLE}.bam -o ${BAM_DIR}/${SAMPLE}_sorted.bam
# Finally index on the sorted aligment file 
samtools index ${BAM_DIR}/${SAMPLE}_sorted.bam

echo "BAM done"

######## lastly, we end with multiqc
conda activate RNA-seq

## Set paths ##
HISAT2_SUMMARY=${SAM_DIR}
MULTIQC_OUT=${SAM_DIR}/multiqc
mkdir -p ${SAM_DIR}

## Run MultiQC ##
echo "Running MultiQC on alignment summary files"

multiqc ${HISAT2_SUMMARY}/*_hisat2_summary.txt -o ${MULTIQC_OUT}

echo "MultiQC complete!"
conda deactivate