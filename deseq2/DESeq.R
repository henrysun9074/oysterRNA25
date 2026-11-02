library(tidyverse)
library(lubridate)
library(viridisLite)
library(readxl)
library(readr)
library(tibble)
library(DESeq2)
library(ggpubr)
library(pheatmap)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(RColorBrewer)
library(clusterProfiler)
library(ashr)
library(vsn)
library(GO.db)
library(ggplot2)

setwd("/work/hs325/bass25")
location_levels <- c("CMT", "DAF", "WC")
time_levels <- c("Jun", "Jul1", "Jul2")

# configure output directory
out_dir <- file.path(getwd(), "deseq2")
figs_dir <- file.path(out_dir, "figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

# read count data
counts_data <- read.table("align/trimmomatic/counts/counts_matrix.txt", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
rownames(counts_data) <- counts_data$Geneid
head(counts_data)

count_matrix <- counts_data[, 7:ncol(counts_data)]

# need to strip the "X.work.hs325.bass25.align.bam." from all column names
colnames(count_matrix) <- gsub("X.work.hs325.bass25.align.bam.", "", colnames(count_matrix))
count_matrix <- as.matrix(count_matrix)
head(count_matrix)

# define trait table using count matrix
colnames(count_matrix)

# example: DAF_YEL3_25Jul25_Gill_S162_L007_sorted
# location = first three letters (factor, DAF | CMT | WC)
# time = middle date (factor, Jun 25-27, Jul 9-11, Jul 23-25)
sample_names <- colnames(count_matrix) %>%
  str_remove("^X\\.work\\.hs325\\.bass25\\.align\\.bam\\.") %>%
  basename() %>%
  str_remove("\\.bam$") %>%
  str_remove("_sorted$")
if (anyDuplicated(sample_names)) stop("Cleaned sample names are duplicated.")
colnames(count_matrix) <- sample_names

# now parse sample metadata, for deseq model design
sample_table <- tibble(sample = sample_names) %>%
  mutate(
    location = str_extract(sample, "^[^_]+"),
    date_string = str_extract(sample, "(?<=_)[0-9]{1,2}(?:Jun|Jul)[0-9]{2}(?=_)"),
    collection_date = dmy(date_string, quiet = TRUE),
    time = case_when(
      between(collection_date, as.Date("2025-06-25"), as.Date("2025-06-27")) ~ "Jun",
      between(collection_date, as.Date("2025-07-09"), as.Date("2025-07-11")) ~ "Jul1",
      between(collection_date, as.Date("2025-07-23"), as.Date("2025-07-25")) ~ "Jul2",
      TRUE ~ NA_character_
    )
  )

# factorize levels of location/date
sample_table <- sample_table %>%
  mutate(
    location = factor(location, levels = location_levels),
    time = factor(time, levels = time_levels),
    group = factor(
      paste(location, time, sep = "_"),
      levels = unlist(lapply(location_levels, function(x) paste(x, time_levels, sep = "_")))
    )
  ) %>%
  dplyr::select(-date_string) %>%
  column_to_rownames("sample")
stopifnot(identical(rownames(sample_table), colnames(count_matrix))) ## check these match

# check should be 3 for all
sample_numbers <- with(sample_table, table(location, time))
print(sample_numbers)

# now save
write.csv(sample_table, file.path(out_dir, "sample_metadata.csv"), row.names = TRUE)

################################################################################

### dds filter - keep only genes with >10 count in >=3 samples

### rerun model if needed
# dds <- DESeqDataSetFromMatrix(
#   countData = count_matrix, colData = sample_table, design = ~ location * time
# )
# dds <- DESeq(dds, betaPrior = FALSE)

# groupsize <- 3L
# keep <- rowSums(counts(dds) >= 10) >= groupsize

# dds_filtered <- DESeqDataSetFromMatrix(
#   countData = counts(dds)[keep, , drop = FALSE],
#   colData = sample_table,
#   design = ~ group
# )
# dds_filtered <- DESeq(dds_filtered, betaPrior = FALSE)
# saveRDS(dds_filtered, file.path(out_dir, "dds_filtered.rds"))

##################### choose here which dds for downstream  ####################

#### both should have already had DESeqDataSetFromMatrix and DESeq ran
######## try again

###### UNFILTERED
# dds_filtered <- readRDS(file.path(out_dir,"dds.rds"))

###### FILTERED
dds_filtered <- readRDS(file.path(out_dir,"dds_filtered.rds"))
stopifnot(identical(design(dds_filtered), ~group)) # make sure this is the filtered ~group fit

################################################################################

####### for plotting
loc_colors <- brewer.pal(n =3, name = "Accent")
time_shapes <- c(8,17,16)
padj.cutoff <- 0.05
lfc.cutoff <- 1 # tested in results() via lfcThreshold, i.e. H0: |LFC| <= lfc.cutoff
top_n <- 100L
gc()

# now do vst for pca plot only
vst_filtered <- varianceStabilizingTransformation(dds_filtered, blind = FALSE)
vst_filtered_matrix <- assay(vst_filtered)
meta <- as.data.frame(colData(dds_filtered))

### interactive 
pca_data <- plotPCA(
  vst_filtered, intgroup = c("location", "time"),
  ntop = min(500L, nrow(vst)), returnData = TRUE
)
percent_var <- round(100 * attr(pca_data, "percentVar"), 2)

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = location, shape = time)) +
  geom_point(size = 4, alpha = 1) +
  scale_color_brewer(palette = "Accent") +
  scale_shape_manual(values = time_shapes) +
  labs(
    x = paste0("PC1: ", percent_var[1], "%"),
    y = paste0("PC2: ", percent_var[2], "%"),
    color = "Location", shape = "Date"
  ) +
  theme_pubr(base_size = 12) +
  theme(text = element_text(face = "bold", family = "sans"),
      legend.text = element_text(face = "plain")) 
print(pca_plot)
ggsave(file.path(figs_dir,"pca.png"), pca_plot, width = 7, height = 7, dpi = 300)

############################### correlation plot
vst_mat <- assay(vst_filtered)  
vst_cor <- cor(vst_mat)  
head(vst_cor)  

# preliminary visualization
heat.colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
hcplot <- pheatmap(vst_cor,
   color = heat.colors,
   fontsize_col = 8,
   fontsize_row = 8) 
hcplot

################ load in model matrix and comparisons to run
model_matrix <- model.matrix(design(dds_filtered), data = meta)
location_mean <- function(loc) {
  colMeans(unique(model_matrix[meta$location == loc, , drop = FALSE]))
}

contrast_list <- list(
  DAF_vs_CMT = location_mean("DAF") - location_mean("CMT"),
  DAF_vs_WC = location_mean("DAF") - location_mean("WC"),
  CMT_vs_WC = location_mean("CMT") - location_mean("WC"),
  DAF_Jun_vs_DAF_Jul2 = c("group", "DAF_Jun", "DAF_Jul2"),
  DAF_Jun_vs_CMT_Jul2 = c("group", "DAF_Jun", "CMT_Jul2"),
  CMT_Jun_vs_CMT_Jul2 = c("group", "CMT_Jun", "CMT_Jul2"),
  DAF_Jun_vs_WC_Jul2 = c("group", "DAF_Jun", "WC_Jul2"),
  CMT_Jun_vs_WC_Jul2 = c("group", "CMT_Jun", "WC_Jul2"),
  WC_Jun_vs_WC_Jul2 = c("group", "WC_Jun", "WC_Jul2")
)

heatmap_samples <- list(
  DAF_vs_CMT = meta$location %in% c("DAF", "CMT"),
  DAF_vs_WC = meta$location %in% c("DAF", "WC"),
  CMT_vs_WC = meta$location %in% c("CMT", "WC"),
  DAF_Jun_vs_DAF_Jul2 = meta$group %in% c("DAF_Jun", "DAF_Jul2"),
  DAF_Jun_vs_CMT_Jul2 = meta$group %in% c("DAF_Jun", "CMT_Jul2"),
  CMT_Jun_vs_CMT_Jul2 = meta$group %in% c("CMT_Jun", "CMT_Jul2"),
  DAF_Jun_vs_WC_Jul2 = meta$group %in% c("DAF_Jun", "WC_Jul2"),
  CMT_Jun_vs_WC_Jul2 = meta$group %in% c("CMT_Jun", "WC_Jul2"),
  WC_Jun_vs_WC_Jul2 = meta$group %in% c("WC_Jun", "WC_Jul2")
)
plot_titles <- c(
  DAF_vs_CMT = "DAF vs CMT",
  DAF_vs_WC = "DAF vs WC",
  CMT_vs_WC = "CMT vs WC",
  DAF_Jun_vs_DAF_Jul2 = "DAF Jun vs DAF Jul2",
  DAF_Jun_vs_CMT_Jul2 = "DAF Jun vs CMT Jul2",
  CMT_Jun_vs_CMT_Jul2 = "CMT Jun vs CMT Jul2",
  DAF_Jun_vs_WC_Jul2 = "DAF Jun vs WC Jul2",
  CMT_Jun_vs_WC_Jul2 = "CMT Jun vs WC Jul2",
  WC_Jun_vs_WC_Jul2 = "WC Jun vs WC Jul2"
)

res_list <- list()
res_shrunk_list <- list()
result_tables <- list()
deg_tables <- list()
volcano_plots <- list()
heatmap_plots <- list()
summary_rows <- list()

############## for each comp - volcano plot, heatmap, DEG list saved
for (comparison in names(contrast_list)) {
  # padj now tests |LFC| > lfc.cutoff, so no post hoc LFC filter is needed
  res <- DESeq2::results(
    dds_filtered, contrast = contrast_list[[comparison]],
    lfcThreshold = lfc.cutoff, altHypothesis = "greaterAbs", alpha = padj.cutoff
  )
  res_list[[comparison]] <- res
  # shrink LFCs so low-count genes don't produce huge, unreliable fold changes;
  # using ashr for shrinkage for combinatorial comparison
  res <- lfcShrink(
    dds_filtered, contrast = contrast_list[[comparison]], res = res, type = "ashr"
  )
  res_shrunk_list[[comparison]] <- res
  tbl <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(
      significance = case_when(
        !is.na(padj) & padj < padj.cutoff & log2FoldChange > 0 ~ "Upregulated",
        !is.na(padj) & padj < padj.cutoff & log2FoldChange < 0 ~ "Downregulated",
        TRUE ~ "Insignificant"
      ),
      significance = factor(significance, levels = c(
        "Downregulated", "Insignificant", "Upregulated"
      ))
    ) %>%
    arrange(padj)
  deg_tbl <- tbl %>% filter(significance != "Insignificant")
  result_tables[[comparison]] <- tbl
  deg_tables[[comparison]] <- deg_tbl
  
  write_csv(deg_tbl, file.path(out_dir, paste0(comparison, "_DEGs.csv")))
  
  summary_rows[[comparison]] <- tibble(
    comparison = comparison,
    n_DEGs = nrow(deg_tbl),
    up_in_first = sum(deg_tbl$log2FoldChange > 0),
    down_in_first = sum(deg_tbl$log2FoldChange < 0)
  )
  
  plot_tbl <- tbl %>%
    filter(!is.na(padj), is.finite(log2FoldChange)) %>%
    mutate(neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)))
  volc <- ggplot(plot_tbl, aes(log2FoldChange, neg_log10_padj, color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-lfc.cutoff, lfc.cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj.cutoff), linetype = "dashed") +
    scale_color_manual(values = c(
      "Downregulated" = "cyan3",
      "Insignificant" = "grey75",
      "Upregulated" = "deeppink4"
      # "Downregulated" = "darkgreen",
      # "Upregulated" = "coral4"
    ), drop = FALSE) +
    theme_pubr(base_size = 12) +
    labs(
      title = unname(plot_titles[comparison]),
      x = "shrunken log2FC",
      y = "-log10-adjusted p-value", color = NULL
    ) +
    theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "plain"),
          text = element_text(face = "bold", family = "sans"),
          legend.text = element_text(face = "plain"))
  volcano_plots[[comparison]] <- volc
  print(volc)
  ggsave(file.path(figs_dir, paste0(comparison, "_volcano.png")), volc,
         width = 8, height = 6, dpi = 300)
  ggsave(file.path(figs_dir, paste0(comparison, "_volcano.pdf")), volc,
         width = 8, height = 6)
  
  # MA plots highlight padj < 0.05 from the |LFC| > lfc.cutoff test (shrunken LFCs shown)
  pdf(file.path(figs_dir, paste0(comparison, "_MA.pdf")), width = 7, height = 6)
  DESeq2::plotMA(res, alpha = padj.cutoff, main = unname(plot_titles[comparison]))
  dev.off()
  
  # heatmap: top significant genes ranked by adjusted p-value
  genes <- head(deg_tbl$gene_id, top_n)
  selected <- which(heatmap_samples[[comparison]])
  selected <- selected[order(meta$location[selected], meta$time[selected])]
  hm_meta <- meta[selected, , drop = FALSE]
  mat <- vst_filtered_matrix[genes, selected, drop = FALSE]
  
  # Remove constant rows before scaling to avoid NA z-scores
  row_sd <- apply(mat, 1, sd)
  mat <- mat[is.finite(row_sd) & row_sd > 0, , drop = FALSE]

  mat <- t(scale(t(mat)))
  annotation <- HeatmapAnnotation(
    Location = hm_meta$location,
    Period = hm_meta$time,
    col = list(
      Location = c("CMT" = "#7FC97F", "DAF" = "#BEAED4", "WC" = "#FDC086"),
      Period = c("Jun" = "lightgoldenrod", "Jul1" = "darkgoldenrod1", "Jul2" = "darkgoldenrod")
    )
  )
  hm <- Heatmap(
    mat, name = "Row z-score",
    col = circlize::colorRamp2(c(-2, 0, 2), rev(brewer.pal(3, "RdBu"))),
    top_annotation = annotation,
    column_title = unname(plot_titles[comparison]),
    show_row_names = FALSE, show_column_names = FALSE,
    cluster_rows = nrow(mat) > 1L, cluster_columns = FALSE
  )
  heatmap_plots[[comparison]] <- hm
  draw(hm)
  pdf(file.path(figs_dir, paste0(comparison, "_top_DEGs_heatmap.pdf")), width = 9, height = 8)
  draw(hm)
  dev.off()
  png(file.path(figs_dir, paste0(comparison, "_top_DEGs_heatmap.png")),
      width = 9, height = 8, units = "in", res = 300)
  draw(hm)
  dev.off()
}

deg_summary <- bind_rows(summary_rows)
print(deg_summary)
write_csv(deg_summary, file.path(out_dir, "DEG_summary.csv"))

# shrunk results are used downstream 
saveRDS(res_shrunk_list, file.path(out_dir, "DEG_contrast_results.rds"))
saveRDS(res_list, file.path(out_dir, "DEG_contrast_results_unshrunk.rds"))

combined_volcano <- ggarrange(
  plotlist = volcano_plots, ncol = 3, nrow = 3,
  common.legend = TRUE, legend = "bottom"
)

print(combined_volcano)

ggsave(file.path(figs_dir, "combined_volcano.pdf"), combined_volcano,
       width = 9, height = 10)
ggsave(file.path(figs_dir, "combined_volcano.png"), combined_volcano,
       width = 9, height = 10)
gc()

# plots for basic QC DESeq object

# run this to get pre shrinkage
# res <- DESeq2::results(
#   dds_filtered, contrast = contrast_list[[comparison]],
#   lfcThreshold = lfc.cutoff, altHypothesis = "greaterAbs", alpha = padj.cutoff
# )

# vsd plot
ntd <- normTransform(dds_filtered)
meanSdPlot(assay(ntd))

# MA plot
plotMA(res, ylim=c(-2,2))

# dispersion plot
plotDispEsts(dds_filtered)

# 

##########################################################

## load in deg summary
deg_df <- read.csv("deseq2/DEG_summary.csv")
deg_results <- readRDS(file.path(out_dir,"DEG_contrast_results.rds"))

##### now check: which genes are significant across multiple comparisons
all_results <- bind_rows(lapply(names(deg_results), function(comp) {
  as.data.frame(deg_results[[comp]]) %>%
    rownames_to_column("gene_id") %>%
    mutate(comparison = comp)
})) %>%
  mutate(
    # padj already tests |LFC| > lfc.cutoff
    significant = !is.na(padj) & padj < padj.cutoff
  )

# count how many comparisons each gene is significant in
shared_genes <- all_results %>%
  filter(significant) %>%
  group_by(gene_id) %>%
  summarise(
    n_significant = n_distinct(comparison),
    .groups = "drop"
  ) %>%
  filter(n_significant >= 2)

# show these genes' results in every comparison
shared_gene_results <- all_results %>%
  inner_join(shared_genes, by = "gene_id") %>%
  dplyr::select(
    gene_id, n_significant, comparison,
    log2FoldChange, pvalue, padj, significant
  ) %>%
  arrange(desc(n_significant), gene_id, comparison)
shared_gene_results

# get only the top genes with greatest abs(LFC) and lowest p-value
# print gene IDs
# output table of comparisons significant in
top_k <- 20L
top_genes <- shared_gene_results %>%
  filter(significant) %>%
  group_by(gene_id) %>%
  summarise(
    n_significant = n_distinct(comparison),
    significant_comparisons = paste(unique(comparison), collapse = "; "),
    max_abs_LFC = max(abs(log2FoldChange)),
    min_pvalue = min(pvalue, na.rm = TRUE),
    min_padj = min(padj),
    .groups = "drop"
  ) %>%
  mutate(
    effect_rank = min_rank(dplyr::desc(max_abs_LFC)),
    pvalue_rank = min_rank(min_padj),
    combined_rank = effect_rank + pvalue_rank
  ) %>%
  arrange(combined_rank, min_padj, dplyr::desc(max_abs_LFC)) %>%
  slice_head(n = top_k)

cat(top_genes$gene_id, sep = "\n")

# Table of significant comparisons and strongest results
top_gene_table <- top_genes %>%
  dplyr::select(
    gene_id, n_significant, significant_comparisons,
    max_abs_LFC, min_pvalue, min_padj
  )
print(top_gene_table, n = Inf, width = Inf)


############################ run GO enrichment ##################################

### GO terms in file /work/hs325/bass25/deseq2/go/GCF_053477285.1-RS_2025_12_gene_ontology.gaf

go_dir <- file.path(out_dir, "go")
gaf_file <- file.path(go_dir, "GCF_053477285.1-RS_2025_12_gene_ontology.gaf")
location_comparisons <- c("DAF_vs_WC", "CMT_vs_WC", "DAF_vs_CMT")
go_comparisons <- unique(c("DAF_Jun_vs_CMT_Jul2", "DAF_Jun_vs_WC_Jul2", location_comparisons))

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

go_results <- list()
go_plots <- list()
for (comparison in go_comparisons) {
  res_tbl <- as.data.frame(deg_results[[comparison]]) %>%
    rownames_to_column("gene_id")
  # background = genes actually tested (independent filtering leaves padj = NA)
  universe <- res_tbl %>% filter(!is.na(padj)) %>% pull(gene_id)
  deg_sig <- res_tbl %>%
    filter(!is.na(padj), padj < padj.cutoff) # padj already tests |LFC| > lfc.cutoff
  gene_sets <- list(
    all = deg_sig$gene_id,
    up = deg_sig$gene_id[deg_sig$log2FoldChange > 0],
    down = deg_sig$gene_id[deg_sig$log2FoldChange < 0]
  )

  for (direction in names(gene_sets)) {
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

      # lollipop of top 15 terms by nominal p (few DEGs rarely survive FDR)
      if (direction == "all") {
        lolli_tbl <- go_results[[key]] %>%
          slice_min(pvalue, n = 15, with_ties = FALSE) %>%
          mutate(
            neg_log10_p = -log10(pvalue),
            Description = str_trunc(Description, 60),
            Description = fct_reorder(Description, neg_log10_p)
          )
        go_plots[[key]] <- ggplot(lolli_tbl, aes(neg_log10_p, Description)) +
          geom_segment(aes(x = 0, xend = neg_log10_p, yend = Description), color = "grey60") +
          geom_point(aes(size = Count, color = p.adjust)) +
          geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
          scale_color_viridis_c(direction = -1, name = "FDR") +
          scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
          theme_pubr(base_size = 12) +
          labs(
            title = paste0(plot_titles[comparison], " - GO ", ont),
            x = "-log10(p-value)", y = NULL, size = "DEGs in term"
          ) +
          theme(plot.title = element_text(size = 12, face = "plain"),
                text = element_text(face = "bold", family = "sans"),
                legend.text = element_text(face = "plain"),
                legend.position = "right")
        print(go_plots[[key]])
        ggsave(file.path(figs_dir, paste0(key, "_GO_lollipop.pdf")), go_plots[[key]],
               width = 8, height = 6)
      }
    }
  }
}

go_table <- bind_rows(go_results) %>% arrange(comparison, direction, ontology, pvalue)
write_csv(go_table, file.path(go_dir, "GO_enrichment_results.csv"))

# quick summary: terms at p < 0.05 and FDR < 0.05 per comparison/direction/ontology
go_table %>%
  group_by(comparison, direction, ontology) %>%
  summarise(n_p05 = sum(pvalue < 0.05), n_fdr05 = sum(p.adjust < 0.05), .groups = "drop") %>%
  print(n = Inf)

############## combined GO lollipop for location-only comparisons
top_terms_per_comp <- 10
loc_go <- go_table %>%
  filter(comparison %in% location_comparisons, direction == "all")
loc_top_ids <- loc_go %>%
  group_by(comparison, ontology) %>%
  slice_min(pvalue, n = top_terms_per_comp, with_ties = FALSE) %>%
  ungroup() %>%
  distinct(ontology, ID)

combined_go_tbl <- loc_go %>%
  semi_join(loc_top_ids, by = c("ontology", "ID")) %>%
  mutate(
    neg_log10_p = -log10(pvalue),
    comparison = factor(comparison, levels = location_comparisons,
                        labels = unname(plot_titles[location_comparisons])),
    ontology = factor(ontology, levels = c("BP", "MF", "CC")),
    Description = str_trunc(Description, 60)
  ) %>%
  group_by(ontology, Description) %>%
  mutate(best_p = max(neg_log10_p)) %>%
  ungroup() %>%
  mutate(Description = fct_reorder(Description, best_p))

combined_go_plot <- ggplot(combined_go_tbl, aes(neg_log10_p, Description)) +
  geom_segment(aes(x = 0, xend = neg_log10_p, yend = Description), color = "grey60") +
  geom_point(aes(size = Count, color = p.adjust)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  facet_grid(ontology ~ comparison, scales = "free_y", space = "free_y") +
  scale_color_viridis_c(direction = -1, name = "FDR") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_pubr(base_size = 11) +
  labs(
    x = "-log10(p-value)", y = NULL, size = "DEGs in term"
  ) +
  theme(plot.title = element_text(size = 12, face = "plain"),
        text = element_text(face = "bold", family = "sans"),
        legend.text = element_text(face = "plain"),
        legend.position = "right",
        strip.background = element_rect(fill = "grey90", color = NA),
        panel.border = element_rect(fill = NA, color = "grey70"))
print(combined_go_plot)

ggsave(file.path(figs_dir, "location_comparisons_GO_lollipop.pdf"), combined_go_plot,
       width = 13, height = 9)
ggsave(file.path(figs_dir, "location_comparisons_GO_lollipop.png"), combined_go_plot,
       width = 13, height = 9, dpi = 300)