#!/bin/bash -e
#SBATCH --job-name=index
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/hs325/bass25/misc/index.out
#SBATCH --error=/work/hs325/bass25/misc/index.err
#SBATCH --partition=scavenger
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

# source conda
source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh

## Load module - HISAT2 already loaded on DCC ## 
module load HISAT2

########## index the thing
GENOME=/work/hs325/cvpan/asms/yu25/ncbi_dataset/data/GCF_053477285.1/GCF_053477285.1_ASM5347728v1_genomic.fna
GTF=/work/hs325/cvpan/asms/yu25/ncbi_dataset/data/GCF_053477285.1/GCF_053477285.1_ASM5347728v1_genomic.gtf
INDEX_DIR=/work/hs325/bass25/align/index
mkdir -p ${INDEX_DIR}

## Create HFM (Hierarchical FM) index ## 
    # Aligns reads to a single reference genome 
hisat2-build \
    -p ${SLURM_CPUS_PER_TASK} \
    ${GENOME} ${INDEX_DIR}/c.virginica_yale25_HFM_index 
