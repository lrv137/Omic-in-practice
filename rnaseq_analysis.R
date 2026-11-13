.libPaths(c(file.path(Sys.getenv("HOME"), "R", "library"), .libPaths()))

# Load required libraries/packages and create output folders

library(DESeq2)
library(data.table)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)


# Record software environments
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")


# Import and validate counts and metadata
counts_df <- fread("data/gene_counts.tsv", header = TRUE, data.table = FALSE, check.names = FALSE)
metadata <- fread("metadata/sample_metadata.tsv", data.table = FALSE)
metadata$batch <- factor(metadata$batch)

gene_ids <- counts_df[[1]]
counts <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(counts) <- gene_ids

mode(counts) <- "numeric"

colnames(counts) <- trimws(as.character(colnames(counts)))
metadata$sample <- trimws(as.character(metadata$sample))


message("--- DEBUG INFO ---")
message("Counts columns:   ", paste(colnames(counts), collapse = ", "))
message("Metadata samples: ", paste(metadata$sample, collapse = ", "))
message("------------------")


#Perform explicit checks before continuing
stopifnot(!anyDuplicated(rownames(counts)))
stopifnot(!anyDuplicated(metadata$sample))
stopifnot(all(counts >= 0))
stopifnot(all(counts == round(counts)))
stopifnot(setequal(colnames(counts), metadata$sample))

metadata <- metadata[match(colnames(counts), metadata$sample), ]
stopifnot(identical(colnames(counts), metadata$sample))
rownames(metadata) <- metadata$sample


# Initial summaries
dim(counts)
summary(colSums(counts))
table(metadata$condition)
head(metadata)


# Plot library sizes before normalization
lib_df <- data.frame(sample = colnames(counts),
    library_size = colSums(counts),
    condition = metadata$condition)

p_lib <- ggplot(lib_df, aes(sample, library_size, fill = condition)) +
    geom_col() +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Total raw counts", title = "Library sizes") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/library_sizes.pdf", p_lib, width = 8, height = 5)


# Create DESeq2 dataset
metadata$condition <- factor(metadata$condition)
metadata$condition <- relevel(metadata$condition, ref = "Control")


# For a simple two-group study
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ condition
)


# If a known batch variable must be accounted for and the study design permits it
# Example only; confirm with the instructor
# metadata$batch <- factor(metadata$batch)
# dds <- DESeqDataSetFromMatrix(counts, metadata, design = ~ batch + condition


# Pre-filter very low-count genes
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
cat("Genes retained:", nrow(dds), "\n")


# Save the prepared object
saveRDS(dds, "results/dds_prefiltered.rds")


# Transform counts for visualization
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)


# Compute and draw PCA
pca_data <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, shape = batch, label = name)) +
    geom_point(size = 3) +
    geom_text_repel(show.legend = FALSE, max.overlaps = Inf) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    theme_bw() +
    labs(title = "PCA of variance-stabilized counts")

ggsave("figures/PCA_samples.pdf", p_pca, width = 7, height = 5.5)
fwrite(pca_data, "results/PCA_coordinates.tsv", sep = "\t")


# Sample-to-sample distance heatmap
sample_dists <- dist(t(vsd_mat))
sample_dist_mat <- as.matrix(sample_dists)
rownames(sample_dist_mat) <- colnames(vsd_mat)
colnames(sample_dist_mat) <- colnames(vsd_mat)

ann_col <- metadata[, c("condition", "batch"), drop = FALSE]

pdf("figures/sample_distance_heatmap.pdf", width = 7, height = 6)
pheatmap(sample_dist_mat, annotation_col = ann_col, annotation_row = ann_col,
    main = "Sample-to-sample distances")
dev.off()


# Read and match the gene set selected for the heatmap
selected_genes <- scan("metadata/selected_genes.txt", what = character(), quiet = TRUE)
present <- intersect(selected_genes, rownames(vsd_mat))
missing <- setdiff(selected_genes, rownames(vsd_mat))
writeLines(missing, "results/selected_genes_not_found.txt")

if (length(present) < 2) stop("Fewer than two selected genes were found")


# Create a z-score heatmap
heat_mat <- vsd_mat[present, , drop = FALSE]
heat_z <- t(scale(t(heat_mat)))
heat_z <- heat_z[complete.cases(heat_z), , drop = FALSE]

ann_col <- metadata[, c("condition", "batch"), drop = FALSE]

pdf("figures/selected_gene_heatmap_zscore.pdf", width = 8, height = 7)
pheatmap(heat_z, annotation_col = ann_col,
    cluster_rows = TRUE, cluster_cols = TRUE,
    show_colnames = TRUE,
    main = "Selected genes: row z-scores")
dev.off()


# Keep samples in experimental order
sample_order <- metadata$sample[order(metadata$condition, metadata$batch)]
heat_z_ordered <- heat_z[, sample_order, drop = FALSE]

pdf("figures/selected_gene_heatmap_ordered.pdf", width = 8, height = 7)
pheatmap(heat_z_ordered, annotation_col = ann_col[sample_order, , drop = FALSE],
    cluster_rows = TRUE, cluster_cols = FALSE,
    main = "Selected genes by experimental group")
dev.off()


# Plot transformed values without row scaling
pdf("figures/selected_gene_heatmap_VST.pdf", width = 8, height = 7)
pheatmap(heat_mat, annotation_col = ann_col, scale = "none",
    main = "Selected genes: VST expression")
dev.off()



# Run DESeq2 model on the raw-count DESeqDataSet
dds <- DESeq(dds)
resultsNames(dds)


# For Treatment compared with Control
res <- results(
    dds,
    contrast = c("condition", "Treatment", "Control"),
    alpha = 0.05
)

summary(res)


# Shrink log2 fold changes for ranking and visualization
res_shrunk <- lfcShrink(
    dds,
    coef = "condition_Treatment_vs_Control",
    res = res,
    type = "apeglm"
)


# Export results
res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

fwrite(res_df, "results/DE_Treatment_vs_Control.tsv", sep = "\t", na = "NA")


# Save dds object with the fitted DESeq2 moddel
saveRDS(dds, "results/dds.rds")


# Create an explicit significance classification
padj_cutoff <- 0.05
lfc_cutoff <- 1

res_df$status <- "Not significant"
res_df$status[!is.na(res_df$padj) & res_df$padj < padj_cutoff &
    res_df$log2FoldChange >= lfc_cutoff] <- "Higher in Treatment"
res_df$status[!is.na(res_df$padj) & res_df$padj < padj_cutoff &
    res_df$log2FoldChange <= -lfc_cutoff] <- "Lower in Treatment"

table(res_df$status)


# MA plot
pdf("figures/MA_Treatment_vs_Control.pdf", width = 7, height = 5)
plotMA(res_shrunk, alpha = 0.05, ylim = c(-5, 5),
    main = "Treatment versus Control")
dev.off()


# Create a volcano plot
plot_df <- res_df
plot_df$minus_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
plot_df$minus_log10_padj[is.na(plot_df$padj)] <- NA_real_

label_pool <- subset(plot_df, status != "Not significant" & !is.na(padj))
label_genes <- head(label_pool[order(label_pool$padj), ], 12)

p_volcano <- ggplot(plot_df, aes(log2FoldChange, minus_log10_padj, color = status)) +
    geom_point(alpha = 0.65, size = 1.4, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = 2) +
    scale_color_manual(values = c(
        "Lower in Treatment" = "#377EB8",
        "Not significant" = "grey70",
        "Higher in Treatment" = "#E41A1C")) +
    geom_text_repel(data = label_genes, aes(label = gene_id),
                size = 3, max.overlaps = Inf, show.legend = FALSE) +
    labs(title = "Treatment versus Control",
            x = "Shrunken log2 fold change",
            y = expression(-log[10](adjusted~p)), color = NULL) +
    theme_bw()

ggsave("figures/volcano_Treatment_vs_Control.pdf", p_volcano,
            width = 7.5, height = 6)
ggsave("figures/volcano_Treatment_vs_Control.png", p_volcano,
            width = 7.5, height = 6, dpi = 300)
