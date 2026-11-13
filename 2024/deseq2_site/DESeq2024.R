library(tidyverse)
library(lubridate)
library(viridisLite)
library(tibble)
library(DESeq2)
library(apeglm)
library(ggpubr)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(vsn)
library(GO.db)
library(ggrepel)
library(limma)

# 2024 vs 2025 location effect
#   - both years are combined and the 2025 WC samples are dropped
#   - design is ~ collection + location instead of ~ group, so that the site
#     effect is a single model coefficient
#   - LFC shrinkage is apeglm instead of ashr

setwd("/work/hs325/bass25/2024")

counts_2024_file <- "/work/hs325/bass25/2024/counts/counts_matrix_exon_2024.txt"
counts_2025_file <- "/work/hs325/bass25/align/trimmomatic/counts/counts_matrix_exon_2Nov26.txt"
gaf_file <- "/work/hs325/bass25/deseq2/go/GCF_053477285.1-RS_2025_12_gene_ontology.gaf"

location_levels <- c("CMT", "DAF") # CMT is the reference, so +LFC = higher in DAF
collection_levels <- c("2024_Jun", "2024_Jul", "2024_Aug", "2025_Jun", "2025_Jul1", "2025_Jul2")

# configure output directory
out_dir <- file.path(getwd(), "deseq2_site")
figs_dir <- file.path(out_dir, "figs")
go_dir <- file.path(out_dir, "go")
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################

# read count data
read_counts <- function(path) {
  x <- read.table(path, header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE,
                  check.names = FALSE, quote = "", comment.char = "")
  rownames(x) <- x$Geneid
  x
}
counts_2024 <- read_counts(counts_2024_file)
counts_2025 <- read_counts(counts_2025_file)

annot_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
stopifnot(identical(counts_2024[, annot_cols], counts_2025[, annot_cols]))

count_matrix <- as.matrix(cbind(counts_2024[, -(1:6)], counts_2025[, -(1:6)]))

# example: DAF_YEL3_25Jul25_Gill_S162_L007
colnames(count_matrix) <- colnames(count_matrix) %>%
  basename() %>%
  str_remove("_sorted\\.bam$")
if (anyDuplicated(colnames(count_matrix))) stop("Cleaned sample names are duplicated.")

# drop WC (2025 only)
count_matrix <- count_matrix[, !str_detect(colnames(count_matrix), "^WC_"), drop = FALSE]
stopifnot(ncol(count_matrix) == 36L, nrow(count_matrix) == nrow(counts_2024))
sample_names <- colnames(count_matrix)

# location = DAF | CMT
# year     = 2024 (SJ) | 2025 (YEL)
# collection = sampling event; CMT and DAF were collected within 1-2 days of
#              each other at each event, so each event holds 3 CMT + 3 DAF
sample_table <- tibble(sample = sample_names) %>%
  mutate(
    location = str_extract(sample, "^[^_]+"),
    date_string = str_extract(sample, "(?<=_)[0-9]{1,2}(?:Jun|Jul|Aug)[0-9]{2}(?=_)"),
    collection_date = dmy(date_string, quiet = TRUE),
    year = factor(year(collection_date)),
    collection = case_when(
      between(collection_date, as.Date("2024-06-20"), as.Date("2024-06-21")) ~ "2024_Jun",
      between(collection_date, as.Date("2024-07-22"), as.Date("2024-07-23")) ~ "2024_Jul",
      between(collection_date, as.Date("2024-08-19"), as.Date("2024-08-20")) ~ "2024_Aug",
      between(collection_date, as.Date("2025-06-25"), as.Date("2025-06-27")) ~ "2025_Jun",
      between(collection_date, as.Date("2025-07-09"), as.Date("2025-07-11")) ~ "2025_Jul1",
      between(collection_date, as.Date("2025-07-23"), as.Date("2025-07-25")) ~ "2025_Jul2",
      TRUE ~ NA_character_
    ),
    location = factor(location, levels = location_levels),
    collection = factor(collection, levels = collection_levels)
  ) %>%
  dplyr::select(-date_string) %>%
  column_to_rownames("sample")
stopifnot(identical(rownames(sample_table), colnames(count_matrix))) ## check these match

# year must line up with the genetic line in the sample name
stopifnot(all(str_detect(sample_names, "_SJ[0-9]_") == (sample_table$year == "2024")))
stopifnot(all(str_detect(sample_names, "_YEL[0-9]_") == (sample_table$year == "2025")))

# check should be 3 for all
sample_numbers <- with(sample_table, table(location, collection))
print(sample_numbers)
stopifnot(all(sample_numbers == 3L), !anyNA(sample_table$collection), !anyNA(sample_table$location))

# now save
write.csv(sample_table, file.path(out_dir, "sample_metadata.csv"), row.names = TRUE)

################################################################################

### dds filter - keep only genes with >=10 counts in >=3 samples
# applied after dropping WC, so genes expressed only in WC are not retained
# groupsize <- 3L
# keep <- rowSums(count_matrix >= 10) >= groupsize
# cat("Genes kept:", sum(keep), "of", nrow(count_matrix), "\n")
# 
# # collection absorbs year/line and within-season differences, leaving
# # location as the DAF vs CMT effect shared across all six events
# dds_filtered <- DESeqDataSetFromMatrix(
#   countData = count_matrix[keep, , drop = FALSE],
#   colData = sample_table,
#   design = ~ collection + location
# )
# dds_filtered <- DESeq(dds_filtered, betaPrior = FALSE)
# saveRDS(dds_filtered, file.path(out_dir, "dds_filtered.rds"))

## load in thing
dds_filtered <- readRDS(file.path(out_dir, "dds_filtered.rds"))

################################################################################

####### for plotting
loc_colors <- c("CMT" = "#7FC97F", "DAF" = "#BEAED4")
collection_shapes <- setNames(c(1, 2, 0, 16, 17, 15), collection_levels) # open = 2024, filled = 2025
padj.cutoff <- 0.05
lfc.cutoff <- 1 
top_n <- 100L
gc()

# now do vst for pca plot only
vst_filtered <- varianceStabilizingTransformation(dds_filtered, blind = FALSE)
vst_filtered_matrix <- assay(vst_filtered)
meta <- as.data.frame(colData(dds_filtered))

pca_data <- plotPCA(
  vst_filtered, intgroup = c("location", "collection"),
  ntop = min(500L, nrow(vst_filtered)), returnData = TRUE
)
percent_var <- round(100 * attr(pca_data, "percentVar"), 2)
stopifnot(identical(pca_data$name, rownames(sample_table)))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = location, shape = collection)) +
  geom_point(size = 4, alpha = 1, stroke = 1) +
  scale_color_manual(values = loc_colors) +
  scale_shape_manual(values = collection_shapes) +
  labs(
    x = paste0("PC1: ", percent_var[1], "%"),
    y = paste0("PC2: ", percent_var[2], "%"),
    color = "Location", shape = "Collection"
  ) +
  theme_pubr(base_size = 12) +
  theme(text = element_text(face = "bold", family = "sans"),
        legend.text = element_text(face = "plain"))
print(pca_plot)
ggsave(file.path(figs_dir, "pca.png"), pca_plot, width = 8, height = 7, dpi = 300,
       device = ragg::agg_png)

# make pca plot with labels
pca_plot_labeled <- pca_plot +
  geom_text_repel(
    aes(label = name),
    size = 2.5, show.legend = FALSE,
    max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0,
    segment.color = "grey50", segment.size = 0.3
  )
print(pca_plot_labeled)
ggsave(file.path(figs_dir, "pca_labeled.png"), pca_plot_labeled, width = 10, height = 9, dpi = 300,
       device = ragg::agg_png)

############################### correlation plot
vst_cor <- cor(vst_filtered_matrix)
heat.colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
hcplot <- pheatmap::pheatmap(vst_cor,
   color = heat.colors,
   fontsize_col = 6,
   fontsize_row = 6,
   silent = TRUE)
ragg::agg_png(file.path(figs_dir, "sample_correlation.png"),
              width = 10, height = 10, units = "in", res = 300)
grid::grid.draw(hcplot$gtable)
dev.off()

gc()
################################################################################

comparison <- "DAF_vs_CMT"
plot_title <- "DAF vs CMT"
site_coef <- "location_DAF_vs_CMT"
stopifnot(site_coef %in% resultsNames(dds_filtered))

# Wald test against LFC = 0
res_unshrunk <- DESeq2::results(
  dds_filtered, name = site_coef, alpha = padj.cutoff
)

res_shrunk <- lfcShrink(
  dds_filtered, coef = site_coef, res = res_unshrunk, type = "apeglm"
)
stopifnot(identical(res_shrunk$padj, res_unshrunk$padj))

result_tbl <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    significance = case_when(
      !is.na(padj) & padj < padj.cutoff & log2FoldChange > lfc.cutoff ~ "Upregulated",
      !is.na(padj) & padj < padj.cutoff & log2FoldChange < -lfc.cutoff ~ "Downregulated",
      TRUE ~ "Insignificant"
    ),
    significance = factor(significance, levels = c(
      "Downregulated", "Insignificant", "Upregulated"
    ))
  ) %>%
  arrange(padj)
deg_tbl <- result_tbl %>% filter(significance != "Insignificant")

write_csv(result_tbl, file.path(out_dir, paste0(comparison, "_all_genes.csv")))
write_csv(deg_tbl, file.path(out_dir, paste0(comparison, "_DEGs.csv")))

deg_summary <- tibble(
  comparison = comparison,
  n_tested = sum(!is.na(res_shrunk$padj)),
  n_DEGs = nrow(deg_tbl),
  up_in_DAF = sum(deg_tbl$log2FoldChange > 0),
  down_in_DAF = sum(deg_tbl$log2FoldChange < 0)
)
print(deg_summary)
write_csv(deg_summary, file.path(out_dir, "DEG_summary.csv"))

# shrunk results are used downstream
saveRDS(res_shrunk, file.path(out_dir, "DEG_contrast_results.rds"))
saveRDS(res_unshrunk, file.path(out_dir, "DEG_contrast_results_unshrunk.rds"))

############## volcano
plot_tbl <- result_tbl %>%
  filter(!is.na(padj), is.finite(log2FoldChange)) %>%
  mutate(neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)))
volc <- ggplot(plot_tbl, aes(log2FoldChange, neg_log10_padj, color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-lfc.cutoff, lfc.cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj.cutoff), linetype = "dashed") +
  # geom_text_repel(
  #   data = ~ filter(.x, significance != "Insignificant"),
  #   aes(label = gene_id),
  #   size = 2.5, show.legend = FALSE,
  #   max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0,
  #   segment.color = "grey50", segment.size = 0.3
  # ) +
  scale_color_manual(values = c(
    "Downregulated" = "#7FC97F",
    "Insignificant" = "grey55",
    "Upregulated" = "#BEAED4"
  ), drop = FALSE) +
  # ylim(0,20)+
  theme_pubr(base_size = 12) +
  labs(
    title = plot_title,
    x = "log2FC",
    y = "-log10-adjusted p-value", color = NULL
  ) +
  theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "plain"),
        text = element_text(face = "bold", family = "sans"),
        legend.text = element_text(face = "plain"))
print(volc)
ggsave(file.path(figs_dir, paste0(comparison, "_volcano.png")), volc,
       width = 8, height = 6, dpi = 300, device = ragg::agg_png)
ggsave(file.path(figs_dir, paste0(comparison, "_volcano.pdf")), volc,
       width = 8, height = 6)

############## heatmap: top significant genes ranked by adjusted p-value
genes <- head(deg_tbl$gene_id, top_n)
selected <- order(meta$location, meta$collection)
hm_meta <- meta[selected, , drop = FALSE]

# remove the collection effect for display only
vst_display <- limma::removeBatchEffect(
  vst_filtered_matrix, batch = meta$collection,
  design = model.matrix(~ location, data = meta)
)
mat <- vst_display[genes, selected, drop = FALSE]

# remove constant rows before scaling to avoid NA z-scores
row_sd <- apply(mat, 1, sd)
mat <- mat[is.finite(row_sd) & row_sd > 0, , drop = FALSE]
mat <- t(scale(t(mat)))

annotation <- HeatmapAnnotation(
  Location = hm_meta$location,
  Year = hm_meta$year,
  Collection = hm_meta$collection,
  col = list(
    Location = loc_colors,
    Year = c("2024" = "grey80", "2025" = "grey30"),
    Collection = setNames(brewer.pal(6, "YlOrBr"), collection_levels)
  )
)
hm <- Heatmap(
  mat, name = "Row z-score",
  col = circlize::colorRamp2(c(-2, 0, 2), rev(brewer.pal(3, "RdBu"))),
  top_annotation = annotation,
  column_title = plot_title,
  show_row_names = FALSE, show_column_names = FALSE,
  cluster_rows = nrow(mat) > 1L, cluster_columns = FALSE
)
draw(hm)
pdf(file.path(figs_dir, paste0(comparison, "_top_DEGs_heatmap.pdf")), width = 9, height = 8)
draw(hm)
dev.off()
ragg::agg_png(file.path(figs_dir, paste0(comparison, "_top_DEGs_heatmap.png")),
              width = 9, height = 8, units = "in", res = 300)
draw(hm)
dev.off()
gc()

############## basic QC for the DESeq object
# MA plot, pre shrinkage
DESeq2::plotMA(res_unshrunk, ylim = c(-2, 2), main = "Unshrunk")
# MA plot, apeglm
DESeq2::plotMA(res_shrunk, ylim = c(-2, 2), main = "apeglm")
# mean-sd plot
ntd <- normTransform(dds_filtered)
print(meanSdPlot(assay(ntd), plot = FALSE)$gg)
# dispersion plot
plotDispEsts(dds_filtered)
median(mcols(dds_filtered)$dispGeneEst)

############################ run GO enrichment ##################################

# parse gaf file
gaf <- read_tsv(gaf_file, comment = "!", col_names = FALSE, col_types = cols(.default = "c"),
                quote = "")
gene2go <- gaf %>%
  transmute(gene_id = X3, qualifier = X4, go_id = X5, ontology = X9) %>%
  filter(!str_detect(qualifier, "^NOT")) %>%
  mutate(ontology = recode(ontology, P = "BP", F = "MF", C = "CC")) %>%
  distinct(gene_id, go_id, ontology)

# GAF only holds the most specific terms, so propagate each annotation to all ancestor terms
ancestor_maps <- list(BP = as.list(GOBPANCESTOR), MF = as.list(GOMFANCESTOR), CC = as.list(GOCCANCESTOR))
gene2go_prop <- bind_rows(lapply(names(ancestor_maps), function(ont) {
  direct <- gene2go %>% filter(ontology == ont)
  anc <- ancestor_maps[[ont]][direct$go_id]
  tibble(
    gene_id = rep(direct$gene_id, lengths(anc) + 1L),
    go_id = unlist(Map(c, direct$go_id, anc), use.names = FALSE),
    ontology = ont
  )
})) %>%
  filter(go_id != "all") %>%
  distinct()

# get GO term names
term_names <- AnnotationDbi::select(GO.db, keys = unique(gene2go_prop$go_id),
                                    columns = "TERM", keytype = "GOID") %>%
  dplyr::rename(go_id = GOID, term = TERM)

# check that the GAF symbols match the count-matrix gene IDs, this should be >0
cat("Genes in dds with >=1 GO term:",
    sum(rownames(dds_filtered) %in% gene2go$gene_id), "of", nrow(dds_filtered), "\n")

# background = genes actually tested (independent filtering leaves padj = NA)
universe <- result_tbl %>% filter(!is.na(padj)) %>% pull(gene_id)
gene_sets <- list(
  all = deg_tbl$gene_id,
  up = deg_tbl$gene_id[deg_tbl$log2FoldChange > 0],
  down = deg_tbl$gene_id[deg_tbl$log2FoldChange < 0]
)

go_results <- list()
for (direction in names(gene_sets)) {
  if (length(gene_sets[[direction]]) == 0L) next
  for (ont in c("BP", "MF", "CC")) {
    t2g <- gene2go_prop %>% filter(ontology == ont) %>% dplyr::select(go_id, gene_id)
    ego <- enricher(
      gene = gene_sets[[direction]], universe = universe,
      TERM2GENE = t2g, TERM2NAME = term_names,
      pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
      minGSSize = 5, maxGSSize = 500
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next
    key <- paste(comparison, direction, ont, sep = "_")
    go_results[[key]] <- as.data.frame(ego) %>%
      mutate(comparison = comparison, direction = direction, ontology = ont, .before = 1)
  }
}
stopifnot(length(go_results) > 0L)

go_table <- bind_rows(go_results) %>% arrange(comparison, direction, ontology, pvalue)
write_csv(go_table, file.path(go_dir, "GO_enrichment_results.csv"))

# quick summary: terms at p < 0.05 and FDR < 0.05 per direction/ontology
go_table %>%
  group_by(comparison, direction, ontology) %>%
  summarise(n_p05 = sum(pvalue < 0.05), n_fdr05 = sum(p.adjust < 0.05), .groups = "drop") %>%
  print(n = Inf)

############## combined GO lollipop by direction, top terms by nominal p
top_terms_per_comp <- 10
direction_labels <- c(all = "All DEGs", up = "Up in DAF", down = "Down in DAF")

dir_top_ids <- go_table %>%
  group_by(direction, ontology) %>%
  slice_min(pvalue, n = top_terms_per_comp, with_ties = FALSE) %>%
  ungroup() %>%
  distinct(ontology, ID)

dir_go_tbl <- go_table %>%
  semi_join(dir_top_ids, by = c("ontology", "ID")) %>%
  mutate(
    neg_log10_p = -log10(pvalue),
    direction = factor(direction, levels = names(direction_labels),
                       labels = unname(direction_labels)),
    ontology = factor(ontology, levels = c("BP", "MF", "CC")),
    Description = str_trunc(Description, 60)
  ) %>%
  group_by(ontology, Description) %>%
  mutate(best_p = max(neg_log10_p)) %>%
  ungroup() %>%
  mutate(Description = fct_reorder(Description, best_p))

direction_go_plot <- ggplot(dir_go_tbl, aes(neg_log10_p, Description)) +
  geom_segment(aes(x = 0, xend = neg_log10_p, yend = Description), color = "grey60") +
  geom_point(aes(size = Count, color = p.adjust)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  facet_grid(ontology ~ direction, scales = "free_y", space = "free_y", drop = FALSE) +
  scale_color_viridis_c(direction = -1, name = "FDR") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_pubr(base_size = 11) +
  labs(
    title = plot_title,
    x = "-log10(p-value)", y = NULL, size = "DEGs in term"
  ) +
  theme(plot.title = element_text(size = 12, face = "plain"),
        text = element_text(face = "bold", family = "sans"),
        legend.text = element_text(face = "plain"),
        legend.position = "right",
        strip.background = element_rect(fill = "grey90", color = NA),
        panel.border = element_rect(fill = NA, color = "grey70"))
print(direction_go_plot)

ggsave(file.path(figs_dir, paste0(comparison, "_GO_direction_lollipop.pdf")),
       direction_go_plot, width = 13, height = 9)
ggsave(file.path(figs_dir, paste0(comparison, "_GO_direction_lollipop.png")),
       direction_go_plot, width = 13, height = 9, dpi = 300, device = ragg::agg_png)

################################################################################

############## 2024 (diploid) vs 2025 (triploid): average counts per gene ##############
# uses the raw, unfiltered count_matrix 
year_vec <- sample_table$year

# per-sample totals
sample_totals <- tibble(
  sample = colnames(count_matrix),
  year = year_vec,
  location = sample_table$location,
  total_counts = colSums(count_matrix),
  genes_detected = colSums(count_matrix > 0),
  mean_counts_per_gene = colMeans(count_matrix)
)
sample_totals %>%
  group_by(year) %>%
  summarise(n = n(),
            mean_total = mean(total_counts), sd_total = sd(total_counts),
            mean_counts_per_gene = mean(mean_counts_per_gene),
            mean_genes_detected = mean(genes_detected), .groups = "drop") %>%
  print()
print(wilcox.test(mean_counts_per_gene ~ year, data = sample_totals))
write_csv(sample_totals, file.path(out_dir, "year_sample_totals.csv"))

# per-gene means within each year both raw and CPM
cpm_matrix <- sweep(count_matrix, 2, colSums(count_matrix), "/") * 1e6
gene_year_tbl <- tibble(
  gene_id = rownames(count_matrix),
  mean_raw_2024 = rowMeans(count_matrix[, year_vec == "2024", drop = FALSE]),
  mean_raw_2025 = rowMeans(count_matrix[, year_vec == "2025", drop = FALSE]),
  mean_cpm_2024 = rowMeans(cpm_matrix[, year_vec == "2024", drop = FALSE]),
  mean_cpm_2025 = rowMeans(cpm_matrix[, year_vec == "2025", drop = FALSE])
) %>%
  mutate(log2_ratio_raw = log2((mean_raw_2025 + 1) / (mean_raw_2024 + 1)),
         log2_ratio_cpm = log2((mean_cpm_2025 + 0.1) / (mean_cpm_2024 + 0.1)))
write_csv(gene_year_tbl, file.path(out_dir, "year_gene_mean_counts.csv"))

# summary over genes expressed in either year
expressed <- gene_year_tbl %>% filter(mean_raw_2024 >= 10 | mean_raw_2025 >= 10)
cat("Genes with mean >= 10 counts in either year:", nrow(expressed), "\n")
expressed %>%
  summarise(mean_raw_2024 = mean(mean_raw_2024), mean_raw_2025 = mean(mean_raw_2025),
            median_raw_2024 = median(mean_raw_2024), median_raw_2025 = median(mean_raw_2025),
            median_log2_ratio_raw = median(log2_ratio_raw),
            median_log2_ratio_cpm = median(log2_ratio_cpm)) %>%
  print(width = Inf)
print(wilcox.test(expressed$mean_raw_2025, expressed$mean_raw_2024, paired = TRUE))

# plot per-sample mean counts per gene and gene-wise 2024 vs 2025 means
year_colors <- c("2024" = "cyan4", "2025" = "cadetblue4")
p_samples <- ggplot(sample_totals, aes(year, mean_counts_per_gene)) +
  geom_boxplot(aes(fill = year), outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(fill = year), width = 0.15, size = 3, shape = 21, color = "white",
              stroke = 0.6, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 24, size = 4,
               fill = "white", color = "black", stroke = 1) +
  ggsignif::geom_signif(comparisons = list(c("2024", "2025")), test = "wilcox.test",
                        map_signif_level = FALSE, textsize = 4, tip_length = 0.01) +
  scale_fill_manual(values = year_colors) +
  labs(x = NULL, y = "Mean raw counts per gene (per sample)") +
  theme_pubr(base_size = 12) + theme(legend.position = "none")
p_samples
