#!/bin/bash
#SBATCH --job-name=read_stats_trimmomatic
#SBATCH --time=06:00:00
#SBATCH --output=/work/hs325/bass25/2024/misc/read_stats_trimmomatic_%a.out
#SBATCH --error=/work/hs325/bass25/2024/misc/read_stats_trimmomatic_%a.err
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## generate summary statistics of the sequencing output for the publication
# usage:  bash read_stats_trimmomatic.sh           (submits fastp array + dependent summary job)
#         sbatch read_stats_trimmomatic.sh summarize   (only re-run the summary)

set -eo pipefail

source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate RNA-seq

## Set paths ##
RAW_DIR=/work/hs325/bass25/2024/raw
TRIMMED_DIR=/work/hs325/bass25/2024/trimmed
ALIGNED_DIR=/work/hs325/bass25/2024/aligned
COUNTS_SUMMARY=/work/hs325/bass25/2024/counts/counts_matrix_exon_2Nov26.txt.summary
FACILITY_METRICS=/work/hs325/bass25/2024/raw/Quality_Metrics.csv
STATS_DIR=/work/hs325/bass25/2024/read_stats
LOG_DIR=/work/hs325/bass25/2024/misc
mkdir -p ${STATS_DIR}/trimmed_fastp ${STATS_DIR}/raw_fastp

# Make list of raw sample names
SAMPLES=($(ls ${RAW_DIR}/*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' | xargs -n 1 basename))
N_SAMPLES=${#SAMPLES[@]}

MODE="${1:-}"

### submit the per-sample array, then the summary once every task has finished
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" && "${MODE}" != "summarize" ]]; then
    script_path="$(readlink -f "${BASH_SOURCE[0]}")"
    array_job="$(sbatch --parsable --array="1-${N_SAMPLES}" "$script_path")"
    array_job="${array_job%%;*}"
    summary_job="$(sbatch --parsable --dependency=afterok:${array_job} \
        --job-name=read_stats_summary --cpus-per-task=1 --mem=4G --time=00:30:00 \
        --output=${LOG_DIR}/read_stats_summary.out --error=${LOG_DIR}/read_stats_summary.err \
        "$script_path" summarize)"
    echo "Submitted fastp array ${array_job} (${N_SAMPLES} samples) and summary job ${summary_job%%;*}"
    exit 0
fi

### per-sample: fastp in report-only mode (no -o/-O, so nothing is written but the reports)
if [[ "${MODE}" != "summarize" ]]; then
    SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}
    echo "Processing sample: " ${SAMPLE}

    fastp -i ${TRIMMED_DIR}/${SAMPLE}_R1.paired.fastq.gz \
          -I ${TRIMMED_DIR}/${SAMPLE}_R2.paired.fastq.gz \
          --disable_adapter_trimming \
          --disable_quality_filtering \
          --disable_length_filtering \
          --disable_trim_poly_g \
          -j ${STATS_DIR}/trimmed_fastp/${SAMPLE}.json \
          -h ${STATS_DIR}/trimmed_fastp/${SAMPLE}.html \
          --thread ${SLURM_CPUS_PER_TASK}

    # raw stats normally come from the poly-G step; only re-scan raw reads if that report is missing
    if [[ ! -s ${TRIMMED_DIR}/${SAMPLE}.polyG.fastp.json ]]; then
        echo "No poly-G fastp report for ${SAMPLE}, scanning raw reads"
        fastp -i ${RAW_DIR}/${SAMPLE}_R1_001.fastq.gz \
              -I ${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz \
              --disable_adapter_trimming \
              --disable_quality_filtering \
              --disable_length_filtering \
              --disable_trim_poly_g \
              -j ${STATS_DIR}/raw_fastp/${SAMPLE}.json \
              -h ${STATS_DIR}/raw_fastp/${SAMPLE}.html \
              --thread ${SLURM_CPUS_PER_TASK}
    fi

    echo "fastp report for " ${SAMPLE} "complete!"
    exit 0
fi

### summary: combine everything into per-sample and overall tables
echo "Summarizing read statistics for ${N_SAMPLES} samples"
python3 - "${RAW_DIR}" "${TRIMMED_DIR}" "${ALIGNED_DIR}" "${COUNTS_SUMMARY}" \
          "${FACILITY_METRICS}" "${STATS_DIR}" "${SAMPLES[@]}" <<'EOF'
import csv, json, os, re, statistics, sys

raw_dir, trim_dir, align_dir, counts_summary, facility_metrics, stats_dir = sys.argv[1:7]
samples = sys.argv[7:]


def fastp_stats(path, prefix, section="before_filtering"):
    d = json.load(open(path))
    s = d["summary"][section]
    # mean of the per-cycle mean quality curves (R1 and R2); approximate for trimmed reads
    # because later cycles hold fewer bases
    curves = [d[f"read{r}_{section}"]["quality_curves"]["mean"] for r in (1, 2)]
    return {
        f"{prefix}_read_pairs": s["total_reads"] // 2,
        f"{prefix}_bases": s["total_bases"],
        f"{prefix}_mean_len_R1": s["read1_mean_length"],
        f"{prefix}_mean_len_R2": s["read2_mean_length"],
        f"{prefix}_mean_cycle_quality": statistics.mean(curves[0] + curves[1]),
        f"{prefix}_Q20_pct": 100 * s["q20_rate"],
        f"{prefix}_Q30_pct": 100 * s["q30_rate"],
        f"{prefix}_GC_pct": 100 * s["gc_content"],
        f"{prefix}_duplication_pct": 100 * d["duplication"]["rate"],
    }


def key_value_file(path):
    out = {}
    for line in open(path):
        if ":" in line:
            k, v = line.split(":", 1)
            out[k.strip()] = float(v)
    return out


def hisat2_stats(path):
    txt = open(path).read()
    pats = {
        "aligned_input_pairs": r"(\d+) reads; of these:",
        "concordant_0_pct": r"\(([\d.]+)%\) aligned concordantly 0 times",
        "concordant_1_pct": r"\(([\d.]+)%\) aligned concordantly exactly 1 time",
        "concordant_multi_pct": r"\(([\d.]+)%\) aligned concordantly >1 times",
        "overall_alignment_pct": r"([\d.]+)% overall alignment rate",
    }
    return {k: float(re.search(p, txt).group(1)) for k, p in pats.items()}


# facility metrics: sum over lanes and both reads, matched by sample-name prefix
facility = {}
if os.path.exists(facility_metrics):
    for row in csv.DictReader(open(facility_metrics)):
        f = facility.setdefault(row["SampleID"], {"yield": 0, "q30": 0, "qsum": 0})
        f["yield"] += int(row["Yield"]); f["q30"] += int(row["YieldQ30"])
        f["qsum"] += int(row["QualityScoreSum"])

# featureCounts: one column per BAM
fc = {}
with open(counts_summary) as fh:
    rows = list(csv.reader(fh, delimiter="\t"))
bams = [os.path.basename(b).replace("_sorted.bam", "") for b in rows[0][1:]]
for i, s in enumerate(bams):
    status = {r[0]: int(r[i + 1]) for r in rows[1:]}
    total = sum(status.values())
    fc[s] = {
        "fc_total_fragments": total,
        "fc_assigned_fragments": status["Assigned"],
        "fc_assigned_pct": 100 * status["Assigned"] / total,
        "fc_multimapping_pct": 100 * status.get("Unassigned_MultiMapping", 0) / total,
        "fc_nofeatures_pct": 100 * status.get("Unassigned_NoFeatures", 0) / total,
        "fc_unmapped_pct": 100 * status.get("Unassigned_Unmapped", 0) / total,
        "fc_singleton_pct": 100 * status.get("Unassigned_Singleton", 0) / total,
    }

per_sample = []
for s in samples:
    row = {"sample": s}

    # raw reads: poly-G step input = raw reads; fall back to raw re-scan
    polyg = os.path.join(trim_dir, f"{s}.polyG.fastp.json")
    if os.path.exists(polyg):
        row.update(fastp_stats(polyg, "raw"))
        after = json.load(open(polyg))["summary"]["after_filtering"]
        row["polyG_bases_removed_pct"] = 100 * (1 - after["total_bases"] / row["raw_bases"])
        # fastp drops reads trimmed to length 0 (all poly-G) even with length filtering off
        row["polyG_dropped_pairs"] = row["raw_read_pairs"] - after["total_reads"] // 2
        row["polyG_dropped_pct"] = 100 * row["polyG_dropped_pairs"] / row["raw_read_pairs"]
    else:
        row.update(fastp_stats(os.path.join(stats_dir, "raw_fastp", f"{s}.json"), "raw"))

    fid = next((k for k in facility if s.startswith(k + "_")), None)
    if fid:
        f = facility[fid]
        row["facility_mean_quality"] = f["qsum"] / f["yield"]
        row["facility_Q30_pct"] = 100 * f["q30"] / f["yield"]

    t = key_value_file(os.path.join(trim_dir, f"{s}.trimmomatic_summary.txt"))
    row.update({
        "trim_input_pairs": int(t["Input Read Pairs"]),
        "trim_both_surviving_pairs": int(t["Both Surviving Reads"]),
        "trim_both_surviving_pct": t["Both Surviving Read Percent"],
        "trim_forward_only_pct": t["Forward Only Surviving Read Percent"],
        "trim_reverse_only_pct": t["Reverse Only Surviving Read Percent"],
        "trim_dropped_pct": t["Dropped Read Percent"],
    })

    row.update(fastp_stats(os.path.join(stats_dir, "trimmed_fastp", f"{s}.json"), "trimmed"))
    row["bases_retained_pct"] = 100 * row["trimmed_bases"] / row["raw_bases"]
    row["total_pairs_retained_pct"] = 100 * row["trimmed_read_pairs"] / row["raw_read_pairs"]

    row.update(hisat2_stats(os.path.join(align_dir, f"{s}_hisat2_summary.txt")))
    row.update(fc.get(s, {}))
    per_sample.append(row)

cols = list(dict.fromkeys(k for r in per_sample for k in r))
with open(os.path.join(stats_dir, "read_stats_per_sample.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(per_sample)

# overall: mean, SD, min, max (and total for count columns)
summary = {}
with open(os.path.join(stats_dir, "read_stats_summary.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric", "n", "mean", "sd", "min", "max", "total"])
    for c in cols[1:]:
        v = [r[c] for r in per_sample if c in r]
        is_count = c.endswith(("_pairs", "_bases", "_fragments"))
        summary[c] = dict(n=len(v), mean=statistics.mean(v), sd=statistics.stdev(v) if len(v) > 1 else 0,
                          min=min(v), max=max(v), total=sum(v) if is_count else None)
        m = summary[c]
        w.writerow([c, m["n"], f"{m['mean']:.4f}", f"{m['sd']:.4f}", f"{m['min']:.4f}",
                    f"{m['max']:.4f}", "" if m["total"] is None else m["total"]])


# short text report with publication-style sentences
def msd(c, scale=1, d=1):
    m = summary[c]
    return (f"{m['mean'] / scale:.{d}f} ± {m['sd'] / scale:.{d}f} "
            f"(range {m['min'] / scale:.{d}f}–{m['max'] / scale:.{d}f})")


lines = [
    f"Samples: {len(per_sample)}",
    f"Raw sequencing: {summary['raw_read_pairs']['total'] / 1e6:.1f} M read pairs in total "
    f"({summary['raw_bases']['total'] / 1e9:.1f} Gb); per sample {msd('raw_read_pairs', 1e6)} M pairs.",
    f"Raw read length: R1 {msd('raw_mean_len_R1', d=0)} bp, R2 {msd('raw_mean_len_R2', d=0)} bp.",
    f"Raw Q30: {msd('raw_Q30_pct')}%; Q20: {msd('raw_Q20_pct')}%; GC: {msd('raw_GC_pct')}%; "
    f"mean per-cycle quality: {msd('raw_mean_cycle_quality')}.",
]
if "facility_mean_quality" in summary:
    lines.append(f"Facility-reported mean quality: {msd('facility_mean_quality', d=2)}; "
                 f"Q30: {msd('facility_Q30_pct')}%.")
lines += [
    f"Duplication (fastp estimate, raw): {msd('raw_duplication_pct')}%.",
    f"Poly-G trimming removed {msd('polyG_bases_removed_pct', d=2)}% of bases and "
    f"{msd('polyG_dropped_pct', d=2)}% of read pairs (entirely poly-G)." if "polyG_bases_removed_pct" in summary else "",
    f"Trimmomatic: {msd('trim_both_surviving_pct', d=2)}% of input read pairs retained as pairs "
    f"(dropped {msd('trim_dropped_pct', d=2)}%).",
    f"Overall: {msd('total_pairs_retained_pct', d=2)}% of raw read pairs retained after all trimming.",
    f"After trimming: {summary['trimmed_read_pairs']['total'] / 1e6:.1f} M read pairs in total; "
    f"per sample {msd('trimmed_read_pairs', 1e6)} M pairs; {msd('bases_retained_pct')}% of bases retained.",
    f"Trimmed read length: R1 {msd('trimmed_mean_len_R1')} bp, R2 {msd('trimmed_mean_len_R2')} bp; "
    f"Q30: {msd('trimmed_Q30_pct')}%.",
    f"HISAT2 overall alignment: {msd('overall_alignment_pct')}%; "
    f"concordant unique: {msd('concordant_1_pct')}%; concordant multi: {msd('concordant_multi_pct')}%.",
]
if "fc_assigned_pct" in summary:
    lines.append(f"featureCounts: {msd('fc_assigned_pct')}% of fragments assigned to genes "
                 f"({msd('fc_assigned_fragments', 1e6)} M per sample); "
                 f"multimapping {msd('fc_multimapping_pct')}%, no feature {msd('fc_nofeatures_pct')}%.")
report = "\n".join(l for l in lines if l)
open(os.path.join(stats_dir, "read_stats_report.txt"), "w").write(report + "\n")
print(report)
EOF

echo "Read statistics written to ${STATS_DIR}"
