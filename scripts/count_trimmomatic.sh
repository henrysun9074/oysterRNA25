#!/bin/bash -e
#SBATCH --job-name=featureCounts_trimmomatic
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/hs325/bass25/2024/misc/featureCounts_trimmomatic.out
#SBATCH --error=/work/hs325/bass25/2024/misc/featureCounts_trimmomatic.err
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## Load module ##
module load Subread

## Set paths ## 
GENOME=/work/hs325/cvpan/asms/yu25/ncbi_dataset/data/GCF_053477285.1
BAM_DIR="/work/hs325/bass25/2024/aligned/bam"
COUNT_DIR="/work/hs325/bass25/2024/counts"
mkdir -p ${COUNT_DIR}

#### fix the genome gtf file
awk 'BEGIN {FS="\t"; OFS="\t"} 
$3 == "exon" && $9 ~ /gene_id ""/ {
    # Extract the transcript_id to use as a backup
    match($9, /transcript_id "[^"]+"/, t);
    if (t[0] != "") {
        sub(/gene_id ""/, "gene_id " substr(t[0], 15), $9);
    } else {
        # If no transcript_id exists, use the product name
        match($9, /product "[^"]+"/, p);
        if (p[0] != "") {
            sub(/gene_id ""/, "gene_id " substr(p[0], 9), $9);
        }
    }
} {print}' ${GENOME}/GCF_053477285.1_ASM5347728v1_genomic.gtf > ${GENOME}/fixed_ncbi_annotation_yale25.gtf

echo "Running featureCounts on all samples simultaneously..."

featureCounts \
    -T ${SLURM_CPUS_PER_TASK} \
    -p --countReadPairs -B \
    -t exon \
    -g gene_id \
    -a ${GENOME}/fixed_ncbi_annotation_yale25.gtf \
    -o ${COUNT_DIR}/counts_matrix_exon_2024.txt \
    ${BAM_DIR}/*_sorted.bam 