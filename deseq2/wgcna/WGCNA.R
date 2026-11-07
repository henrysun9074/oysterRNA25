# WGCNA on dds_filtered.rds
# adapted from Callie script
################################################################################
# BiocManager::install("DESeq2", force = TRUE)  # one-time; do not rerun every session
# BiocManager::install("GenomicRanges", force = TRUE)
# BiocManager::install("preprocessCore", force = TRUE, update = FALSE)
# BiocManager::install("clusterProfiler", force = TRUE, update = FALSE)

library(tidyverse)
library(DESeq2)
library(WGCNA)
library(matrixStats)
library(gridExtra)
library(patchwork)
library(ggpubr)
library(clusterProfiler)
library(GO.db)

options(stringsAsFactors = FALSE)
allowWGCNAThreads() 

setwd("/work/hs325/bass25")
location_levels <- c("CMT", "DAF", "WC")
time_levels <- c("Jun", "Jul1", "Jul2")

out_dir <- file.path(getwd(), "deseq2")
wgcna_dir <- file.path(out_dir, "wgcna")
wgcna_figs <- file.path(wgcna_dir, "figs")
dir.create(wgcna_figs, recursive = TRUE, showWarnings = FALSE)

#### parameters ####
gene_fraction <- 0.1   # top fraction of genes by VST variance (1 = all filtered genes)
soft_power <- NULL      # NULL = auto pick lowest power with signed R^2 >= r2_cutoff
r2_cutoff <- 0.80
min_module_size <- 50
merge_cut_height <- 0.25 # merge modules with eigengene correlation > 0.75
run_tag <- paste0("top", gene_fraction, "var")

# Load data
dds_filtered <- readRDS(file.path(out_dir, "dds_filtered.rds"))
stopifnot(identical(design(dds_filtered), ~group))

# load metadata straight from the dds object
coldata <- as.data.frame(colData(dds_filtered)) %>%
  mutate(
    location = factor(location, levels = location_levels),
    time = factor(time, levels = time_levels)
  )

# Normalize (VST) and select genes
dds_norm <- vst(dds_filtered, blind = TRUE) # blind = TRUE for unsupervised analysis
norm.counts <- assay(dds_norm)

gene_variance <- rowVars(norm.counts)
names(gene_variance) <- rownames(norm.counts)
number_of_wgcna_genes <- ceiling(length(gene_variance) * gene_fraction)
cat("Genes used for WGCNA:", number_of_wgcna_genes, "of", length(gene_variance), "\n")

selected_wgcna_genes <- names(sort(gene_variance, decreasing = TRUE))[seq_len(number_of_wgcna_genes)]

# # WGCNA wants samples in rows, genes in columns
datExpr <- t(norm.counts[selected_wgcna_genes, , drop = FALSE])
storage.mode(datExpr) <- "numeric"
stopifnot(identical(rownames(datExpr), rownames(coldata)))
write.csv(datExpr, file.path(wgcna_dir, paste0(run_tag, "_norm.counts.csv")))

# # QC: outlier genes/samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  cat("Removing", sum(!gsg$goodGenes), "genes and", sum(!gsg$goodSamples), "samples\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  coldata <- coldata[gsg$goodSamples, , drop = FALSE]
}

# sample clustering to eyeball outliers
sample_tree <- hclust(dist(datExpr), method = "average")
pdf(file.path(wgcna_figs, paste0(run_tag, "_sample_clustering.pdf")), width = 12, height = 6)
plot(sample_tree, main = "Sample clustering (VST, WGCNA genes)", sub = "", xlab = "", cex = 0.6)
dev.off()

# Soft-threshold power - automatically calculate
power <- c(1:10, seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(
  data = datExpr,
  powerVector = power,
  networkType = "signed",
  verbose = 5
)
sft.data <- data.frame(sft$fitIndices)
# for signed networks the WGCNA tutorials use -sign(slope) * R^2
sft.data$signed.R.sq <- -sign(sft.data$slope) * sft.data$SFT.R.sq

R.sq.plot <- ggplot(sft.data, aes(Power, signed.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.05) +
  geom_hline(yintercept = r2_cutoff, color = "red") +
  labs(x = "Power", y = "SFT.R.sq") +
  theme_pubr()

mean.k.plot <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.05 * max(sft.data$mean.k.)) +
  labs(x = "Power", y = "Mean connectivity") +
  theme_pubr()

sft_plot <- grid.arrange(R.sq.plot, mean.k.plot, nrow = 2)
ggsave(device = ragg::agg_png, file.path(wgcna_figs, paste0(run_tag, "_soft_threshold.png")), sft_plot, width = 6, height = 8)

if (is.null(soft_power)) {
  ok <- sft.data$Power[sft.data$signed.R.sq >= r2_cutoff]
  if (length(ok) > 0) {
    soft_power <- min(ok)
  } else {
    # WGCNA default for signed networks with < 30 samples
    soft_power <- 18
    warning("No power reached R^2 >= ", r2_cutoff)
  }
}
cat("Soft threshold power:", soft_power, "\n")

# Build network + modules
temp_cor <- cor
cor <- WGCNA::cor

bwnetwk <- blockwiseModules(
  datExpr = datExpr,
  maxBlockSize = max(5000, ncol(datExpr)),
  networkType = "signed",
  TOMType = "signed",
  power = soft_power,
  minModuleSize = min_module_size,
  mergeCutHeight = merge_cut_height,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  randomSeed = 1234,
  verbose = 3
)

cor <- temp_cor
saveRDS(bwnetwk, file.path(wgcna_dir, paste0(run_tag, "_wgcna_results.rds")))

############################################
## load existing network
# bwnetwk <- readRDS(file.path(wgcna_dir, paste0(run_tag, "_wgcna_results.rds")))

############################################ downstream analysis
# modules
module_colors <- as.character(unlist(bwnetwk$colors, use.names = FALSE))
stopifnot(length(module_colors) == ncol(datExpr), !anyNA(module_colors))
names(module_colors) <- colnames(datExpr)

module_size_table <- tibble(module = module_colors) %>%
  dplyr::count(module, name = "number_of_genes") %>%
  arrange(desc(number_of_genes))
print(module_size_table)

# dendrogram, unmerged vs merged colors 
pdf(file.path(wgcna_figs, paste0(run_tag, "_module_dendrogram.pdf")), width = 10, height = 6)
for (b in seq_along(bwnetwk$dendrograms)) {
  idx <- bwnetwk$blockGenes[[b]]
  plotDendroAndColors(
    bwnetwk$dendrograms[[b]],
    cbind(bwnetwk$unmergedColors[idx], bwnetwk$colors[idx]),
    c("unmerged", "merged"),
    main = paste("Gene dendrogram, block", b),
    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05
  )
}
dev.off()

# gene -> module table
gene_module_table <- tibble(gene_id = names(module_colors), module = module_colors)
write_csv(gene_module_table, file.path(wgcna_dir, paste0(run_tag, "_gene_module_assignments.csv")))

# get module eigengenes
module_eigengenes <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
module_eigengenes <- orderMEs(module_eigengenes)
rownames(module_eigengenes) <- rownames(datExpr)
stopifnot(identical(rownames(module_eigengenes), rownames(coldata)))
eigengene_columns <- colnames(module_eigengenes)

sample_level_eigengenes <- coldata %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, location, time, group) %>%
  bind_cols(as_tibble(module_eigengenes))
write_csv(sample_level_eigengenes,
          file.path(wgcna_dir, paste0(run_tag, "_sample_level_eigengenes.csv")))

# collapse to the 9 collection events (location x time)
event_level_eigengenes <- sample_level_eigengenes %>%
  group_by(group, location, time) %>%
  summarise(across(all_of(eigengene_columns), mean), RNAseq_samples = n(), .groups = "drop") %>%
  arrange(location, time)
stopifnot(nrow(event_level_eigengenes) == length(location_levels) * length(time_levels))
write_csv(event_level_eigengenes,
          file.path(wgcna_dir, paste0(run_tag, "_event_level_eigengenes.csv")))

# eigengene boxplots by group
me_long <- sample_level_eigengenes %>%
  pivot_longer(all_of(eigengene_columns), names_to = "module", values_to = "eigengene")
me_box <- ggplot(me_long, aes(time, eigengene, fill = location, color = location)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.7,
               alpha = 0.5, color = "grey40") +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 1.8) +
  facet_wrap(~module, scales = "free_y") +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  labs(x = NULL, y = "Module eigengene", fill = "Location", color = "Location") +
  theme_pubr(base_size = 11)
print(me_box)
ggsave(file.path(wgcna_figs, paste0(run_tag, "_eigengene_boxplots.pdf")), me_box,
       width = 12, height = 8)
ggsave(device = ragg::agg_png, file.path(wgcna_figs, paste0(run_tag, "_eigengene_boxplots.png")), me_box,
       width = 12, height = 8)

# Trait matrix with location + time
trait_matrix <- coldata %>%
  transmute(
    location_CMT = as.integer(location == "CMT"),
    location_DAF = as.integer(location == "DAF"),
    location_WC  = as.integer(location == "WC"),
    time_Jun  = as.integer(time == "Jun"),
    time_Jul1 = as.integer(time == "Jul1"),
    time_Jul2 = as.integer(time == "Jul2"),
    time_numeric = as.integer(time) # ordinal Jun < Jul1 < Jul2 trend
  ) %>%
  as.matrix()
rownames(trait_matrix) <- rownames(coldata)

eigengene_matrix <- as.matrix(module_eigengenes)
stopifnot(identical(rownames(eigengene_matrix), rownames(trait_matrix)))

# Module-trait correlations
pairwise_correlation_results <- function(eigengene_matrix, trait_matrix) {
  expand_grid(module = colnames(eigengene_matrix), trait = colnames(trait_matrix)) %>%
    pmap_dfr(function(module, trait) {
      x <- eigengene_matrix[, module]
      y <- trait_matrix[, trait]
      keep <- complete.cases(x, y)
      if (sum(keep) >= 3 && sd(x[keep]) > 0 && sd(y[keep]) > 0) {
        tst <- cor.test(x[keep], y[keep], method = "pearson")
        r <- unname(tst$estimate); p <- tst$p.value
      } else {
        r <- NA_real_; p <- NA_real_
      }
      tibble(module = module, trait = trait, correlation = r, n = sum(keep), p_value = p)
    }) %>%
    mutate(FDR = p.adjust(p_value, method = "BH"),
           absolute_correlation = abs(correlation))
}

module_trait_results <- pairwise_correlation_results(eigengene_matrix, trait_matrix) %>%
  arrange(FDR, p_value, desc(absolute_correlation))
write_csv(module_trait_results,
          file.path(wgcna_dir, paste0(run_tag, "_module_trait_correlations.csv")))
print(module_trait_results, n = 25)

# heatmap of module-trait correlations + module sizes
module_order <- colnames(eigengene_matrix)
trait_order <- colnames(trait_matrix)
trait_display_names <- c(
  location_CMT = "CMT", location_DAF = "DAF", location_WC = "WC",
  time_Jun = "Jun", time_Jul1 = "Jul1", time_Jul2 = "Jul2",
  time_numeric = "Time trend\n(Jun→Jul2)"
)

heatmap_plot_data <- module_trait_results %>%
  mutate(
    module = factor(module, levels = rev(module_order)),
    trait = factor(trait, levels = trait_order),
    significant = !is.na(p_value) & p_value < 0.05,
    p_display = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "<0.001",
      TRUE ~ formatC(p_value, format = "f", digits = 3)
    ),
    cell_label = if_else(
      significant,
      paste0(formatC(correlation, format = "f", digits = 2), "\n(p=", p_display, ")"),
      ""
    )
  )

correlation_heatmap <- ggplot(heatmap_plot_data, aes(trait, module, fill = correlation)) +
  geom_tile(color = "grey85", linewidth = 0.5) +
  geom_text(aes(label = cell_label), size = 3, fontface = "bold", lineheight = 0.9) +
  geom_vline(xintercept = c(3.5, 6.5), color = "grey40") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D73027",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  scale_x_discrete(labels = trait_display_names, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_pubr(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(face = "bold"),
        legend.position = "left",
        plot.margin = margin(5, 5, 5, 5))

module_count_plot_data <- module_size_table %>%
  mutate(module_ME = factor(paste0("ME", module), levels = rev(module_order))) %>%
  filter(!is.na(module_ME))

gene_count_plot <- ggplot(module_count_plot_data, aes(number_of_genes, module_ME)) +
  geom_col(aes(fill = module), width = 0.72, show.legend = FALSE) +
  scale_fill_identity() +
  geom_text(aes(label = scales::comma(number_of_genes)), hjust = -0.15, size = 4) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(limits = c(0, max(module_count_plot_data$number_of_genes) * 1.4),
                     expand = expansion(mult = c(0, 0))) +
  labs(subtitle = "Number of genes", x = NULL, y = NULL) +
  theme_pubr(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = margin(5, 10, 5, 0))

final_wgcna_plot <- correlation_heatmap + gene_count_plot + plot_layout(widths = c(5.5, 1.4))
print(final_wgcna_plot)

ggsave(file.path(wgcna_figs, paste0(run_tag, "_module_trait_heatmap.pdf")), final_wgcna_plot,
       width = 10, height = max(4, 0.45 * length(module_order) + 2), dpi = 600, bg = "white")
ggsave(device = ragg::agg_png, file.path(wgcna_figs, paste0(run_tag, "_module_trait_heatmap.png")), final_wgcna_plot,
       width = 10, height = max(4, 0.45 * length(module_order) + 2))

# Module membership (kME) + gene significance -> hub genes
nSamples <- nrow(datExpr)

# module membership: correlation of each gene with each eigengene (genes x modules)
module_membership <- WGCNA::cor(datExpr, eigengene_matrix, use = "p")
module_membership_p <- corPvalueStudent(module_membership, nSamples)
colnames(module_membership) <- sub("^ME", "kME_", colnames(module_membership))
colnames(module_membership_p) <- sub("^ME", "p.kME_", colnames(module_membership_p))

# gene significance: correlation of each gene with each trait (genes x traits)
gene_significance <- WGCNA::cor(datExpr, trait_matrix, use = "p")
gene_significance_p <- corPvalueStudent(gene_significance, nSamples)
colnames(gene_significance) <- paste0("GS_", colnames(gene_significance))
colnames(gene_significance_p) <- paste0("p.GS_", colnames(gene_significance_p))

gene_info <- tibble(gene_id = colnames(datExpr), module = module_colors) %>%
  bind_cols(as_tibble(module_membership), as_tibble(module_membership_p),
            as_tibble(gene_significance), as_tibble(gene_significance_p))
write_csv(gene_info, file.path(wgcna_dir, paste0(run_tag, "_gene_MM_GS.csv")))

# hub genes: top 20 by kME within each (non-grey) module
hub_genes <- gene_info %>%
  filter(module != "grey") %>%
  rowwise() %>%
  mutate(kME_own = c_across(starts_with("kME_"))[match(paste0("kME_", module),
                                                       colnames(module_membership))]) %>%
  ungroup() %>%
  group_by(module) %>%
  slice_max(kME_own, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(module, gene_id, kME_own, starts_with("GS_"))
write_csv(hub_genes, file.path(wgcna_dir, paste0(run_tag, "_hub_genes_top20.csv")))

# MM vs GS scatter for the module-trait pairs with nominal p < 0.05
sig_pairs <- module_trait_results %>% filter(!is.na(p_value), p_value < 0.05, module != "MEgrey")
pdf(file.path(wgcna_figs, paste0(run_tag, "_MM_vs_GS.pdf")), width = 6, height = 6)
for (i in seq_len(nrow(sig_pairs))) {
  mod <- sub("^ME", "", sig_pairs$module[i])
  trt <- sig_pairs$trait[i]
  in_mod <- module_colors == mod
  verboseScatterplot(
    abs(module_membership[in_mod, paste0("kME_", mod)]),
    abs(gene_significance[in_mod, paste0("GS_", trt)]),
    xlab = paste("Module membership in", mod),
    ylab = paste("Gene significance for", trt),
    main = paste(mod, "module vs", trt, "\n"),
    col = mod, cex.main = 1, cex.lab = 1, cex.axis = 1
  )
}
dev.off()


############################ GO enrichment for selected modules ##################
# same GAF parsing / ancestor propagation / enricher settings as in DEG GO analysis

go_modules <- c("turquoise", "brown", "blue", "green", "black", "red", "pink")
stopifnot(all(go_modules %in% module_colors))

go_dir <- file.path(out_dir, "go")
gaf_file <- file.path(go_dir, "GCF_053477285.1-RS_2025_12_gene_ontology.gaf")

gaf <- read_tsv(gaf_file, comment = "!", col_names = FALSE, col_types = cols(.default = "c"),
                quote = "")
gene2go <- gaf %>%
  transmute(gene_id = X3, qualifier = X4, go_id = X5, ontology = X9) %>%
  filter(!str_detect(qualifier, "^NOT")) %>%
  mutate(ontology = recode(ontology, P = "BP", F = "MF", C = "CC")) %>%
  distinct(gene_id, go_id, ontology)

# propagate GAF annotations to all ancestor terms
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

term_names <- AnnotationDbi::select(GO.db, keys = unique(gene2go_prop$go_id),
                                    columns = "TERM", keytype = "GOID") %>%
  dplyr::rename(go_id = GOID, term = TERM)

# background = all genes that went into the WGCNA NOT entire transcriptome
wgcna_universe <- colnames(datExpr)
cat("WGCNA genes with >=1 GO term:",
    sum(wgcna_universe %in% gene2go$gene_id), "of", length(wgcna_universe), "\n")

module_go_results <- list()
for (mod in go_modules) {
  mod_genes <- names(module_colors)[module_colors == mod]
  for (ont in c("BP", "MF", "CC")) {
    t2g <- gene2go_prop %>% filter(ontology == ont) %>% dplyr::select(go_id, gene_id)
    ego <- enricher(
      gene = mod_genes, universe = wgcna_universe,
      TERM2GENE = t2g, TERM2NAME = term_names,
      pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
      minGSSize = 5, maxGSSize = 500
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next
    module_go_results[[paste(mod, ont, sep = "_")]] <- as.data.frame(ego) %>%
      mutate(module = paste0("ME", mod), ontology = ont, .before = 1)
  }
}

module_go_table <- bind_rows(module_go_results) %>% arrange(module, ontology, pvalue)
write_csv(module_go_table,
          file.path(wgcna_dir, paste0(run_tag, "_module_GO_enrichment.csv")))

module_go_table %>%
  group_by(module, ontology) %>%
  summarise(n_p05 = sum(pvalue < 0.05), n_fdr05 = sum(p.adjust < 0.05), .groups = "drop") %>%
  print(n = Inf)

# one lollipop figure per module, facet strips filled with the module color
top_terms_per_module <- 10
module_go_plots <- list()
for (mod in go_modules) {
  me <- paste0("ME", mod)
  mod_go_tbl <- module_go_table %>%
    filter(module == me) %>%
    group_by(ontology) %>%
    slice_min(pvalue, n = top_terms_per_module, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      neg_log10_p = -log10(pvalue),
      ontology = factor(ontology, levels = c("BP", "MF", "CC")),
      Description = str_trunc(Description, 60),
      Description = fct_reorder(Description, neg_log10_p)
    )
  if (nrow(mod_go_tbl) == 0) next

  # white strip text on dark module colors, black on light ones
  strip_rgb <- col2rgb(mod)[, 1]
  strip_luminance <- sum(strip_rgb * c(0.299, 0.587, 0.114)) / 255
  strip_text_color <- if (strip_luminance < 0.5) "white" else "black"

  module_go_plots[[mod]] <- ggplot(mod_go_tbl, aes(neg_log10_p, Description)) +
    geom_segment(aes(x = 0, xend = neg_log10_p, yend = Description), color = "grey60") +
    geom_point(aes(size = Count, color = p.adjust)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    facet_grid(ontology ~ module, scales = "free_y", space = "free_y") +
    scale_color_viridis_c(direction = -1, name = "FDR") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_pubr(base_size = 11) +
    labs(x = "-log10(p-value)", y = NULL, size = "Module genes\nin term") +
    theme(text = element_text(face = "bold", family = "sans"),
          legend.text = element_text(face = "plain"),
          legend.position = "right",
          strip.background = element_rect(fill = mod, color = NA),
          strip.text = element_text(color = strip_text_color, face = "bold"),
          panel.border = element_rect(fill = NA, color = "grey70"))
  print(module_go_plots[[mod]])

  ggsave(file.path(wgcna_figs, paste0(run_tag, "_", me, "_GO_lollipop.pdf")),
         module_go_plots[[mod]], width = 8, height = 9)
  ggsave(device = ragg::agg_png, file.path(wgcna_figs, paste0(run_tag, "_", me, "_GO_lollipop.png")),
         module_go_plots[[mod]], width = 8, height = 9, dpi = 300)
}
