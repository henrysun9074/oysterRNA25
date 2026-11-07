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

# R/4.4.3 has no WGCNA; R440 = module R/4.4.0 + ~/R/4.4.0-rhel9 (see ~/bin/R440)
cd /work/hs325/bass25/deseq2/wgcna

$HOME/bin/R440 WGCNA.R