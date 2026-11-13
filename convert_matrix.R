#!/usr/bin/env Rscript

# 1. Set up a user-writable personal library path for cluster environments
user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

# 2. Ensure BiocManager and org.Hs.eg.db are installed in your personal library
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = user_lib)
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db", ask = FALSE, lib = user_lib)
}

library(org.Hs.eg.db, lib.loc = user_lib)

input_file <- "data/gene_counts_ensembl.tsv"
output_file <- "data/gene_counts.tsv"

cat("1. Loading count matrix from:", input_file, "\n")
# check.names = FALSE prevents R from prepending 'X' to numeric column names
counts <- read.delim(input_file, row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)

cat("2. Stripping Ensembl ID version numbers...\n")
clean_ids <- sub("\\..*", "", rownames(counts))

cat("3. Mapping Ensembl IDs to Gene Symbols...\n")
symbols <- mapIds(org.Hs.eg.db,
                  keys = clean_ids,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")

# Store symbols as a column
counts$gene_symbol <- ifelse(is.na(symbols), rownames(counts), symbols)

cat("4. Aggregating counts while preserving original row order...\n")
# Track original row positions
counts$original_order <- 1:nrow(counts)

# Capture the index of first appearance for each gene symbol
min_order <- aggregate(original_order ~ gene_symbol, data = counts, FUN = min)

# Aggregate counts (excluding the tracking column from summation)
counts_data_only <- counts[, !colnames(counts) %in% "original_order"]
counts_aggregated <- aggregate(. ~ gene_symbol, data = counts_data_only, FUN = sum)

# Merge back with min_order and re-sort to restore original sequence
counts_aggregated <- merge(counts_aggregated, min_order, by = "gene_symbol")
counts_aggregated <- counts_aggregated[order(counts_aggregated$original_order), ]

# Set row names and clean up helper columns
rownames(counts_aggregated) <- counts_aggregated$gene_symbol
counts_aggregated$gene_symbol <- NULL
counts_aggregated$original_order <- NULL
counts <- counts_aggregated

cat("5. Saving result to:", output_file, "\n")
write.table(counts, output_file, sep = "\t", quote = FALSE, col.names = NA)

cat("Successfully converted and saved with original columns and order preserved!\n")
