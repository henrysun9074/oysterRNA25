#!/bin/bash -e
#SBATCH --job-name=fastp
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/hs325/bass25/misc/fastp.out
#SBATCH --error=/work/hs325/bass25/misc/fastp.err
#SBATCH --partition=scavenger
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate RNA-seq

cd /work/hs325/bass25

for r1 in raw/*_R1_001.fastq.gz; do
    # Define the corresponding R2 file path
    r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    # Extract the base sample name (strips path and suffix)
    base=$(basename "$r1" _R1_001.fastq.gz)
    
    echo "Processing sample: $base"
    
    # Run fastp for paired-end data
    fastp -i "$r1" -I "$r2" \
          -o "clean/${base}_R1.trimmed.fastq.gz" \
          -O "clean/${base}_R2.trimmed.fastq.gz" \
          --detect_adapter_for_pe \
          -h "clean/${base}.fastp.html" \
          -j "clean/${base}.fastp.json" \
          --thread 4
done
