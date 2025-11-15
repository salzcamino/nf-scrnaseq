#!/usr/bin/env Rscript

# Cell Quality Control and Filtering using Seurat
# This script performs QC filtering on single-cell RNA-seq data

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: seurat_qc.R <matrix_dir> <sample_id> [min_genes] [min_cells] [max_mito_pct] [min_counts] [max_counts] [max_genes]")
}

# Parse arguments
matrix_dir <- args[1]
sample_id <- args[2]
min_genes <- ifelse(length(args) >= 3, as.numeric(args[3]), 200)
min_cells <- ifelse(length(args) >= 4, as.numeric(args[4]), 3)
max_mito_pct <- ifelse(length(args) >= 5, as.numeric(args[5]), 5)
min_counts <- ifelse(length(args) >= 6, as.numeric(args[6]), 1000)
max_counts <- ifelse(length(args) >= 7 && args[7] != "null", as.numeric(args[7]), Inf)
max_genes <- ifelse(length(args) >= 8 && args[8] != "null", as.numeric(args[8]), Inf)

cat("Quality Control Parameters:\n")
cat(sprintf("  Min genes per cell: %d\n", min_genes))
cat(sprintf("  Min cells per gene: %d\n", min_cells))
cat(sprintf("  Max mitochondrial %%: %.1f\n", max_mito_pct))
cat(sprintf("  Min UMI counts: %d\n", min_counts))
cat(sprintf("  Max UMI counts: %s\n", ifelse(is.infinite(max_counts), "None", max_counts)))
cat(sprintf("  Max genes: %s\n", ifelse(is.infinite(max_genes), "None", max_genes)))

# Load data
cat("\nLoading data from:", matrix_dir, "\n")
counts <- Read10X(data.dir = matrix_dir)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = min_cells,
    min.features = min_genes
)

cat(sprintf("Initial cells: %d\n", ncol(seurat_obj)))
cat(sprintf("Initial genes: %d\n", nrow(seurat_obj)))

# Calculate mitochondrial percentage
# Try different mitochondrial gene patterns
mito_pattern <- "^MT-"
if (sum(grepl("^mt-", rownames(seurat_obj))) > 0) {
    mito_pattern <- "^mt-"
}

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mito_pattern)

# Calculate ribosomal percentage
ribo_pattern <- "^RPL|^RPS"
if (sum(grepl("^Rpl|^Rps", rownames(seurat_obj))) > 0) {
    ribo_pattern <- "^Rpl|^Rps"
}
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = ribo_pattern)

# Generate pre-filtering plots
cat("\nGenerating pre-filtering QC plots...\n")

p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste0(sample_id, "_qc_pre_filter.pdf"), width = 14, height = 8)
print(p1)
grid.arrange(p2, p3, ncol = 2)
dev.off()

# Apply filters
cat("\nApplying filters...\n")
seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA >= min_genes &
             nFeature_RNA <= max_genes &
             nCount_RNA >= min_counts &
             nCount_RNA <= max_counts &
             percent.mt < max_mito_pct
)

cat(sprintf("Cells after filtering: %d\n", ncol(seurat_obj)))
cat(sprintf("Genes after filtering: %d\n", nrow(seurat_obj)))

# Calculate filtering statistics
filtering_stats <- data.frame(
    metric = c("initial_cells", "initial_genes", "filtered_cells", "filtered_genes",
               "cells_removed", "genes_removed", "percent_cells_kept", "percent_genes_kept"),
    value = c(
        ncol(counts), nrow(counts),
        ncol(seurat_obj), nrow(seurat_obj),
        ncol(counts) - ncol(seurat_obj),
        nrow(counts) - nrow(seurat_obj),
        round(ncol(seurat_obj) / ncol(counts) * 100, 2),
        round(nrow(seurat_obj) / nrow(counts) * 100, 2)
    )
)

write.csv(filtering_stats, paste0(sample_id, "_filtering_stats.csv"), row.names = FALSE)

# Generate post-filtering plots
cat("\nGenerating post-filtering QC plots...\n")

p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste0(sample_id, "_qc_post_filter.pdf"), width = 14, height = 8)
print(p1)
grid.arrange(p2, p3, ncol = 2)
dev.off()

# Save filtered object
cat("\nSaving filtered Seurat object...\n")
saveRDS(seurat_obj, file = paste0(sample_id, "_filtered.rds"))

cat("\nQC filtering complete!\n")
