# Load libraries
library(DESeq2)
library(ggplot2)


# Load the fitted DESeq2 dataset
dds <- readRDS("results/dds.rds")


# Select a gene to investigate (replace GENE_ID with an identifier present in the matrix)
gene_to_plot <- "GENE_ID"
count_df <- plotCounts(dds, gene = gene_to_plot,
                intgroup = c("condition", "batch"), returnData = TRUE)

p_count <- ggplot(count_df, aes(condition, count, color = batch)) +
    geom_point(position = position_jitter(width = 0.08), size = 3) +
    scale_y_log10() + theme_bw() +
    labs(title = gene_to_plot, y = "Normalized count (log scale)")

ggsave(paste0("figures/counts_", gene_to_plot, ".pdf"),
        p_count, width = 5, height = 4.5)
