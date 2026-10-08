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
counts_data <- read.table("align/counts/counts_matrix.txt", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
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
  select(-date_string) %>%
  column_to_rownames("sample")
stopifnot(identical(rownames(sample_table), colnames(count_matrix))) ## check these match

# check should be 3 for all
sample_numbers <- with(sample_table, table(location, time))
print(sample_numbers)

# now save
write.csv(sample_table, file.path(out_dir, "sample_metadata.csv"), row.names = TRUE)

################################################################################

####### for plotting
loc_colors <- brewer.pal(n =3, name = "Accent")
time_shapes <- c(8,17,16)

## build dds object
# dds <- DESeqDataSetFromMatrix(
#   countData = count_matrix, colData = sample_table, design = ~ location * time
# )
# dds <- DESeq(dds, betaPrior = FALSE)

# save and reload
# saveRDS(dds, file = file.path(out_dir,"dds.rds"))
dds <- readRDS(file.path(out_dir,"dds.rds"))

# now do vst for pca plot or load in unfiltered
vst <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_matrix <- assay(vst)
# saveRDS(vst, file = file.path(out_dir,"vst.rds"))
# vst <- readRDS(file.path(out_dir,"vst.rds"))

### interactive 
pca_data <- plotPCA(
  vsd, intgroup = c("location", "time"),
  ntop = min(500L, nrow(vsd)), returnData = TRUE
)
percent_var <- round(100 * attr(pca_data, "percentVar"), 2)

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = location, shape = time)) +
  geom_point(size = 4, alpha = 0.8) +
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

### location 
pca_data_l <- plotPCA(
  vsd, intgroup = "location",
  ntop = min(500L, nrow(vsd)), returnData = TRUE
)
percent_var <- round(100 * attr(pca_data, "percentVar"), 2)

pca_plot_l <- ggplot(pca_data_l, aes(PC1, PC2, color = location)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_brewer(palette = "Accent") +
  labs(
    x = paste0("PC1: ", percent_var[1], "%"),
    y = paste0("PC2: ", percent_var[2], "%"),
    color = "Location"
  ) +
  theme_pubr(base_size = 12) +
  theme(text = element_text(face = "bold", family = "sans"),
        legend.text = element_text(face = "plain")) 
print(pca_plot_l)
ggsave(file.path(figs_dir,"pca_location.png"), pca_plot_l, width = 7, height = 7, dpi = 300)

### time
pca_data_t <- plotPCA(
  vsd, intgroup = "time",
  ntop = min(500L, nrow(vsd)), returnData = TRUE
)
percent_var <- round(100 * attr(pca_data, "percentVar"), 2)

pca_plot_t <- ggplot(pca_data_t, aes(PC1, PC2, shape = time)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_shape_manual(values = time_shapes) +
  labs(
    x = paste0("PC1: ", percent_var[1], "%"),
    y = paste0("PC2: ", percent_var[2], "%"),
    shape = "Date"
  ) +
  theme_pubr(base_size = 12) +
  theme(text = element_text(face = "bold", family = "sans"),
        legend.text = element_text(face = "plain")) 
print(pca_plot_t)
ggsave(file.path(figs_dir,"pca_time.png"), pca_plot_t, width = 7, height = 7, dpi = 300)

#### correlation plot
vst_mat <- assay(vst)  
vst_cor <- cor(vst_mat)  
head(vst_cor)  

## TODO: modify this plot to remove sample names but add in labels for location/date
## Add in legend and save
heat.colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
hcplot <- pheatmap(vst_cor,
   color = heat.colors,
   fontsize_col = 8,
   fontsize_row = 8) 
hcplot
# ggsave(file.path(figs_dir,"clustering_raw.png"), plot = hcplot, width = 7, height = 7, dpi = 300)


################################################################################

### filter before DEG analysis any genes with counts < 10 in >=3 samples
nrow(dds) #35738
groupsize <- 3
keep <- rowSums(counts(dds) >= 10) >= groupsize
head(keep)
dds_filtered <- dds[keep, ]
nrow(dds_filtered) #24136

# save the filtered
saveRDS(dds_filtered, file = file.path(out_dir,"dds_filtered.rds"))

################################################################################

## TODO: fix model to do within location, within a timepoint... rather than global

### start here for DEG analysis
dds <- readRDS(file.path(out_dir,"dds_filtered.rds"))
dds_filtered <- readRDS(file.path(out_dir,"dds_filtered.rds"))

dds <- DESeq(dds)
dds_filtered <- DESeq(dds_filtered)

### use the filtered for now cuz callie said so :p
results <- results(dds_filtered)
summary(results)
results0.05 <- results(dds_filtered, alpha = 0.05)
results0.05 <- results0.05[order(results0.05$padj), ]
summary(results0.05)

# start with more lenient LFC cutoff
padj.cutoff <- 0.05
lfc.cutoff <- 1

nrow(as.data.frame(results0.05))
# only genes with padj <0.05
nrow(filter(as.data.frame(results0.05), padj<padj.cutoff))
# only genes with padj <0.05 and LFC > 1
nrow(filter(as.data.frame(results0.05),padj<padj.cutoff,abs(log2FoldChange)>lfc.cutoff))

### plot dispersion
plotDispEsts(dds_filtered, 
             main = "Dispersion estimates for model design date + location") 

### plot MA (LFC vs normalized counts)
plotMA(results)

### make volcano plot
results_df <- as.data.frame(results) %>% 
  rownames_to_column("gene_id") %>% 
  mutate(
    significance = ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1,
                          "Significant", "Not significant")
  )

### Volcano plot
ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed") +
  scale_color_manual(values = c("red", "blue")) + 
  theme_pubr() +
  labs(
    x = "log2 fold change",
    y = "-log10 adjusted p-value") +
  theme(legend.position = "none")


