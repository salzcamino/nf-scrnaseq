# nf-scrnaseq

A Nextflow pipeline for single-cell RNA-seq analysis using Scanpy.

## Overview

This pipeline performs quality control and analysis of single-cell RNA-sequencing data. Currently implements:

- **Data Import**: Support for multiple input formats (10X Genomics, H5AD, CSV)
- **Quality Control**: Cell and gene filtering based on QC metrics
- **QC Visualization**: Comprehensive plots and reports

## Quick Start

### Prerequisites

- Nextflow >= 22.10.0
- Docker, Singularity, or Conda

### Installation

```bash
# Clone the repository
git clone https://github.com/salzcamino/nf-scrnaseq.git
cd nf-scrnaseq

# Build the Docker image (required for Docker profile)
docker build -t nf-scrnaseq:latest .

# OR, use Conda profile instead (no Docker build needed)
# See Profiles section below

# Run with test data
nextflow run main.nf -profile test,docker
```

### Basic Usage

```bash
# Run with 10X Genomics data
nextflow run main.nf \
  --input /path/to/10x/folder \
  -profile docker

# Run with H5AD file
nextflow run main.nf \
  --input /path/to/data.h5ad \
  --input_format h5ad \
  -profile docker

# Run with custom QC parameters
nextflow run main.nf \
  --input /path/to/data \
  --min_genes 300 \
  --max_genes 5000 \
  --max_pct_mt 10 \
  -profile docker
```

## Input Formats

The pipeline supports the following input formats:

1. **10X Genomics**: Directory containing `matrix.mtx`, `genes.tsv`, `barcodes.tsv`
2. **H5AD**: AnnData HDF5 format (`.h5ad`)
3. **10X H5**: 10X HDF5 format (`.h5`)
4. **CSV/TSV**: Gene expression matrix (genes as rows, cells as columns)

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to input data (directory or file) |

### QC Filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_genes` | 200 | Minimum genes per cell |
| `--min_cells` | 3 | Minimum cells per gene |
| `--max_genes` | 2500 | Maximum genes per cell |
| `--max_counts` | null | Maximum counts per cell (auto if null) |
| `--max_pct_mt` | 5 | Maximum mitochondrial percentage |
| `--exclude_mt` | false | Exclude mitochondrial genes |
| `--exclude_ribo` | false | Exclude ribosomal genes |

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | ./results | Output directory |
| `--publish_dir_mode` | copy | Publishing mode (copy/symlink/move) |

## Output Structure

```
results/
├── import/
│   ├── raw_data.h5ad              # Raw imported data
│   └── import_summary.txt         # Import statistics
└── qc/
    ├── qc_filtered.h5ad          # QC-filtered data
    ├── qc_metrics.csv            # QC metrics per cell
    ├── qc_plots.pdf              # QC visualization plots
    └── qc_summary.txt            # QC filtering summary
```

## Profiles

The pipeline supports multiple execution profiles:

- **`docker`**: Use Docker containers (requires building image first: `docker build -t nf-scrnaseq:latest .`)
- **`conda`**: Use Conda environments (recommended if Docker not available)
- **`singularity`**: Use Singularity containers
- **`test`**: Run with test dataset (combine with docker/conda, e.g., `-profile test,conda`)
  - Uses synthetic test data with 50 cells and 100 genes
  - Automatically applies relaxed QC thresholds (min_genes=10, max_genes=100, min_cells=1)
  - Test data is intentionally sparse to enable quick validation

Examples:
```bash
# Using Docker (after building the image)
nextflow run main.nf --input data/ -profile docker

# Using Conda (no build required)
nextflow run main.nf --input data/ -profile conda

# Test with Conda profile (uses adjusted QC parameters automatically)
nextflow run main.nf -profile test,conda
```

**Note**: For real scRNA-seq data, use the default QC parameters (min_genes=200, max_genes=2500) which are appropriate for typical datasets.

## QC Metrics

The pipeline calculates and visualizes the following QC metrics:

- **n_genes_by_counts**: Number of genes detected per cell
- **total_counts**: Total UMI counts per cell
- **pct_counts_mt**: Percentage of mitochondrial gene counts
- **pct_counts_ribo**: Percentage of ribosomal gene counts

## Visualization

The QC module generates comprehensive plots including:

- Distribution histograms for counts, genes, and MT percentage
- Violin plots for key metrics
- Scatter plots showing relationships between metrics
- Threshold lines indicating filtering cutoffs

## Pipeline Status

**Current Version**: 0.1.0

**Implemented**:
- Data import (10X, H5AD, CSV formats)
- Quality control and filtering
- QC visualization

**Coming Soon**:
- Normalization and scaling
- Highly variable gene selection
- Dimensionality reduction (PCA, UMAP, t-SNE)
- Clustering
- Cell type annotation
- Differential expression analysis

## Requirements

### Software Dependencies

The pipeline uses the following key packages:

- scanpy >= 1.9.3
- anndata >= 0.9.2
- pandas >= 2.0.3
- matplotlib >= 3.7.2
- seaborn >= 0.12.2

### Container/Environment Options

Dependencies can be satisfied in multiple ways:

1. **Docker**: Build the container using the provided `Dockerfile`
   ```bash
   docker build -t nf-scrnaseq:latest .
   ```

2. **Conda**: Use the provided `environment.yml`
   ```bash
   nextflow run main.nf -profile conda ...
   ```
   (Nextflow will automatically create the conda environment)

3. **Singularity**: Convert the Docker image to Singularity format

## Troubleshooting

### WSL (Windows Subsystem for Linux) Issues

If you're running on WSL and encounter Python errors like "Could not find platform independent libraries" or "ModuleNotFoundError: No module named 'io'":

**Solution**: The pipeline is already configured to use `$HOME/.nextflow-conda-cache` for conda environments, which is on the Linux filesystem. If you still have issues:

1. Make sure you're running the pipeline from a Linux directory (not `/mnt/c/...`)
2. Or, manually set the conda cache directory:
   ```bash
   export NXF_CONDA_CACHEDIR=$HOME/.nextflow-conda-cache
   nextflow run main.nf -profile test,conda
   ```

3. Install mamba for faster and more reliable conda environment creation:
   ```bash
   conda install -n base -c conda-forge mamba
   ```

### Conda Environment Creation is Slow

The first run will take longer as conda needs to download and install all dependencies. Subsequent runs will use the cached environment. Using mamba (recommended) significantly speeds this up.

## Help

For detailed help message:

```bash
nextflow run main.nf --help
```

## Citation

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- **Scanpy**: Wolf, F.A., et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology, 19(1), 15.

## License

MIT License

## Contact

For issues and questions, please use the GitHub issue tracker.