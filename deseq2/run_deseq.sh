#!/bin/bash
#SBATCH --job-name=deseq_2025
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=6:00:00
#SBATCH --output=/work/hs325/bass25/deseq2/run_deseq.%j.log

# R440 = module R/4.4.0 + ~/R/4.4.0-rhel9 (skips the container-built ~/R/4.4.0)
cd /work/hs325/bass25/deseq2
$HOME/bin/R440 /work/hs325/bass25/deseq2/DESeq.R
