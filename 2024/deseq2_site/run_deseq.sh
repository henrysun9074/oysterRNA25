#!/bin/bash
#SBATCH --job-name=deseq_site
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00
#SBATCH --output=/work/hs325/bass25/2024/deseq2_site/run_deseq.%j.log

module purge
module load R/4.4.3

Rscript /work/hs325/bass25/2024/deseq2_site/DESeq2024.R