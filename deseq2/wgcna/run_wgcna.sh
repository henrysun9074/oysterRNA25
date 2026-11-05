#!/bin/bash -e
#SBATCH --job-name=wgcna
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/hs325/bass25/misc/wgcna.out
#SBATCH --error=/work/hs325/bass25/misc/wgcna.err
#SBATCH --partition=scavenger
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

module load R/4.4.3

cd /work/hs325/bass25/deseq2/wgcna

Rscript WGCNA.R