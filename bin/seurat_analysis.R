#!/usr/bin/env Rscript

# Seurat Analysis Pipeline
# Normalization, scaling, dimensionality reduction, clustering

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(gridExtra)
    library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: seurat_analysis.R <input_rds> <sample_id> [n_hvgs] [n_pcs] [resolution] [norm_method]")
}

# Parse arguments
input_rds <- args[1]
sample_id <- args[2]
n_hvgs <- ifelse(length(args) >= 3, as.numeric(args[3]), 2000)
n_pcs <- ifelse(length(args) >= 4, as.numeric(args[4]), 50)
resolution <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.8)
norm_method <- ifelse(length(args) >= 6, args[6], "lognorm")

cat("Analysis Parameters:\n")
cat(sprintf("  Number of HVGs: %d\n", n_hvgs))
cat(sprintf("  Number of PCs: %d\n", n_pcs))
cat(sprintf("  Clustering resolution: %.2f\n", resolution))
cat(sprintf("  Normalization method: %s\n", norm_method))

# Load data
cat("\nLoading Seurat object...\n")
seurat_obj <- readRDS(input_rds)

cat(sprintf("Cells: %d\n", ncol(seurat_obj)))
cat(sprintf("Genes: %d\n", nrow(seurat_obj)))

# Normalization
cat("\nNormalizing data...\n")
if (norm_method == "sctransform") {
    seurat_obj <- SCTransform(
        seurat_obj,
        vars.to.regress = "percent.mt",
        verbose = FALSE,
        variable.features.n = n_hvgs
    )
} else {
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = n_hvgs)

    # Scale data
    all_genes <- rownames(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj, features = all_genes, vars.to.regress = "percent.mt")
}

# Plot variable features
cat("\nPlotting highly variable genes...\n")
p_hvg <- VariableFeaturePlot(seurat_obj)
p_hvg_labeled <- LabelPoints(plot = p_hvg, points = head(VariableFeatures(seurat_obj), 10), repel = TRUE)

pdf(paste0(sample_id, "_hvg.pdf"), width = 10, height = 6)
print(p_hvg_labeled)
dev.off()

# PCA
cat("\nRunning PCA...\n")
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = n_pcs)

# Visualize PCA
pdf(paste0(sample_id, "_pca.pdf"), width = 12, height = 8)
print(DimPlot(seurat_obj, reduction = "pca"))
print(DimHeatmap(seurat_obj, dims = 1:9, cells = 500, balanced = TRUE))
print(ElbowPlot(seurat_obj, ndims = n_pcs))
dev.off()

# Determine dimensionality
# Use first n_pcs or fewer based on elbow
use_pcs <- min(n_pcs, 30)  # Default to first 30 PCs

# Clustering
cat(sprintf("\nClustering with resolution %.2f...\n", resolution))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:use_pcs)
seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

cat(sprintf("Number of clusters: %d\n", length(unique(Idents(seurat_obj)))))

# UMAP
cat("\nRunning UMAP...\n")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:use_pcs)

# Visualizations
cat("\nGenerating visualizations...\n")

# UMAP plots
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident", pt.size = 0.5)

pdf(paste0(sample_id, "_umap.pdf"), width = 14, height = 6)
print(p1 + p2)
dev.off()

# Feature plots for key QC metrics
p_qc <- FeaturePlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2)

pdf(paste0(sample_id, "_umap_qc.pdf"), width = 15, height = 5)
print(p_qc)
dev.off()

# Find markers for all clusters
cat("\nFinding cluster markers...\n")
all_markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

write.csv(all_markers, paste0(sample_id, "_markers.csv"), row.names = FALSE)

# Get top markers per cluster
top_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

write.csv(top_markers, paste0(sample_id, "_top_markers.csv"), row.names = FALSE)

# Heatmap of top markers
if (nrow(top_markers) > 0) {
    top_genes <- unique(top_markers$gene)
    if (length(top_genes) > 0) {
        pdf(paste0(sample_id, "_marker_heatmap.pdf"), width = 12, height = 10)
        print(DoHeatmap(seurat_obj, features = top_genes) + NoLegend())
        dev.off()
    }
}

# Generate summary statistics
cluster_stats <- data.frame(
    cluster = names(table(Idents(seurat_obj))),
    n_cells = as.vector(table(Idents(seurat_obj)))
)

write.csv(cluster_stats, paste0(sample_id, "_cluster_stats.csv"), row.names = FALSE)

# Save processed object
cat("\nSaving processed Seurat object...\n")
saveRDS(seurat_obj, file = paste0(sample_id, "_processed.rds"))

# Export to h5ad for compatibility with scanpy
if (requireNamespace("SeuratDisk", quietly = TRUE)) {
    cat("\nExporting to h5ad format...\n")
    library(SeuratDisk)
    SaveH5Seurat(seurat_obj, filename = paste0(sample_id, ".h5seurat"), overwrite = TRUE)
    Convert(paste0(sample_id, ".h5seurat"), dest = "h5ad", overwrite = TRUE)
}

# Generate HTML report summary
cat("\nGenerating summary report...\n")

summary_text <- sprintf("
# Single-cell RNA-seq Analysis Summary
## Sample: %s

### Data Overview
- Total cells: %d
- Total genes: %d
- Number of clusters: %d

### Analysis Parameters
- Highly variable genes: %d
- Principal components: %d
- Clustering resolution: %.2f
- Normalization method: %s

### QC Metrics
- Median genes per cell: %.0f
- Median UMI counts per cell: %.0f
- Median mitochondrial percentage: %.2f%%

### Clustering Results
%s

Analysis completed successfully!
",
    sample_id,
    ncol(seurat_obj),
    nrow(seurat_obj),
    length(unique(Idents(seurat_obj))),
    n_hvgs,
    use_pcs,
    resolution,
    norm_method,
    median(seurat_obj$nFeature_RNA),
    median(seurat_obj$nCount_RNA),
    median(seurat_obj$percent.mt),
    paste(capture.output(print(cluster_stats)), collapse = "\n")
)

writeLines(summary_text, paste0(sample_id, "_summary.txt"))

cat("\nAnalysis complete!\n")
