# Usage Guide

## Table of Contents

1. [Getting Started](#getting-started)
2. [Input Preparation](#input-preparation)
3. [Running the Pipeline](#running-the-pipeline)
4. [Understanding Outputs](#understanding-outputs)
5. [Parameter Selection](#parameter-selection)
6. [Troubleshooting](#troubleshooting)

## Getting Started

### System Requirements

- **Nextflow**: Version 22.10.0 or higher
- **Container Engine**: Docker, Singularity, or Podman (recommended)
  - OR Conda/Mamba (alternative)
- **Hardware**: Minimum 32GB RAM (64GB+ recommended for large datasets)
- **Storage**: ~50GB for STAR index + output space

### Installation

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

## Input Preparation

### Creating a Samplesheet

The samplesheet is a CSV file that tells the pipeline where to find your data:

```csv
sample,fastq_1,fastq_2,protocol
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,10x_3prime
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,10x_3prime
```

**Column descriptions:**
- `sample`: Unique identifier for each sample
- `fastq_1`: Path to read 1 FASTQ (or cDNA reads for 10x)
- `fastq_2`: Path to read 2 FASTQ (cell barcode + UMI for 10x)
- `protocol`: Optional, overrides global `--protocol` parameter

### Reference Genome Preparation

#### Option 1: Pre-built References

For human and mouse, you can use pre-built indices (if available):

```bash
--genome GRCh38  # Human (hg38)
--genome GRCm39  # Mouse (mm39)
```

#### Option 2: Custom Reference

Download FASTA and GTF files:

```bash
# Example: Human genome from GENCODE
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz

gunzip *.gz

# Use in pipeline
--fasta GRCh38.primary_assembly.genome.fa \
--gtf gencode.v43.primary_assembly.annotation.gtf
```

## Running the Pipeline

### Basic Command

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --genome GRCh38 \
  -profile docker
```

### With Custom Parameters

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --fasta genome.fa \
  --gtf genes.gtf \
  --protocol 10x_3prime \
  --min_genes 500 \
  --max_mito_pct 10 \
  --resolution 0.5 \
  --n_hvgs 3000 \
  -profile docker \
  -resume
```

### Using Different Profiles

```bash
# Docker (recommended)
-profile docker

# Singularity (for HPC)
-profile singularity

# Conda
-profile conda

# Test with small dataset
-profile test,docker

# Multiple profiles
-profile docker,test
```

### Running on HPC

```bash
# SLURM
nextflow run main.nf \
  -profile singularity \
  -process.executor slurm \
  -process.queue your_queue

# PBS
nextflow run main.nf \
  -profile singularity \
  -process.executor pbs
```

### Cloud Execution

```bash
# AWS
nextflow run main.nf \
  -profile docker \
  -bucket-dir s3://my-bucket/work \
  -process.executor awsbatch \
  -process.queue my-queue

# Google Cloud
nextflow run main.nf \
  -profile docker \
  -bucket-dir gs://my-bucket/work
```

## Understanding Outputs

### Quality Control

**File**: `results/qc/sample_qc_pre_filter.pdf`
- Violin plots showing nFeature_RNA, nCount_RNA, percent.mt
- Scatter plots of UMI counts vs. genes and mitochondrial percentage
- **Use to**: Determine appropriate QC thresholds

**File**: `results/qc/sample_filtering_stats.csv`
- Number of cells before/after filtering
- Percentage of cells retained
- **Use to**: Assess filtering stringency

### Analysis Results

**File**: `results/analysis/sample_umap.pdf`
- UMAP visualization colored by cluster
- **Use to**: Visualize cell populations and clustering quality

**File**: `results/analysis/sample_markers.csv`
- All marker genes for each cluster
- Columns: cluster, gene, p_val, avg_log2FC, pct.1, pct.2
- **Use to**: Identify cell types

**File**: `results/analysis/sample_processed.rds`
- Complete Seurat object with all analysis results
- **Use to**: Further downstream analysis in R

### MultiQC Report

**File**: `results/multiqc/multiqc_report.html`
- Aggregated QC metrics across all samples
- **Use to**: Compare samples and identify outliers

## Parameter Selection

### Quality Control Thresholds

#### Standard 10x v3 Chemistry (PBMC)
```bash
--min_genes 200 \
--max_mito_pct 5 \
--min_counts 1000
```

#### Strict Quality (High-quality cells)
```bash
--min_genes 500 \
--max_mito_pct 5 \
--min_counts 2000 \
--max_counts 50000 \
--max_genes 6000
```

#### Permissive (Include more cells)
```bash
--min_genes 100 \
--max_mito_pct 10 \
--min_counts 500
```

### Analysis Parameters

#### Standard Analysis
```bash
--n_hvgs 2000 \
--n_pcs 50 \
--resolution 0.8
```

#### High Resolution (More clusters)
```bash
--n_hvgs 3000 \
--n_pcs 50 \
--resolution 1.2
```

#### Low Resolution (Fewer, broader clusters)
```bash
--n_hvgs 2000 \
--n_pcs 30 \
--resolution 0.4
```

### Normalization Methods

#### Log-normalization (Default, faster)
```bash
--normalization_method lognorm
```

#### SCTransform (Better for heterogeneous datasets)
```bash
--normalization_method sctransform
```

## Troubleshooting

### Pipeline Fails to Start

**Error**: "Command not found: nextflow"
```bash
# Ensure Nextflow is installed and in PATH
export PATH=$PATH:/path/to/nextflow
```

**Error**: "Unable to find image"
```bash
# Ensure Docker/Singularity is running
docker --version
singularity --version
```

### Out of Memory Errors

**Solution 1**: Increase max memory
```bash
--max_memory 256.GB
```

**Solution 2**: Process samples sequentially
```bash
# Edit nextflow.config
process.maxForks = 1
```

### STAR Alignment Issues

**Error**: "Genome index does not exist"
```bash
# Ensure genome files are accessible
ls -lh /path/to/genome.fa
ls -lh /path/to/genes.gtf
```

**Error**: "Not enough memory to generate genome"
```bash
# STAR requires ~30GB for human genome
--max_memory 64.GB
```

### No Cells Pass Filtering

**Cause**: QC thresholds too strict

**Solution**: Check pre-filtering QC plots and adjust:
```bash
# Lower thresholds
--min_genes 100 \
--max_mito_pct 10 \
--min_counts 500
```

### Empty Clusters

**Cause**: Resolution too high

**Solution**: Decrease clustering resolution:
```bash
--resolution 0.3
```

### Pipeline Hangs

**Solution**: Check resource availability
```bash
# Monitor system resources
top
htop

# Check Nextflow processes
nextflow log
```

### Resume Failed Run

```bash
# Resume from last successful step
nextflow run main.nf -resume

# Clean work directory if resume fails
rm -rf work/
nextflow run main.nf
```

## Advanced Usage

### Custom Configuration

Create a custom config file `custom.config`:

```nextflow
process {
    withName: STAR_ALIGN {
        cpus = 24
        memory = 128.GB
    }
}
```

Run with custom config:
```bash
nextflow run main.nf -c custom.config
```

### Multiple Samples with Different Parameters

Not directly supported. Process samples separately:

```bash
# Sample 1
nextflow run main.nf --input sample1.csv --outdir results_s1 --min_genes 200

# Sample 2
nextflow run main.nf --input sample2.csv --outdir results_s2 --min_genes 500
```

### Integration with Downstream Tools

#### Export to Scanpy (Python)

```python
import scanpy as sc

# Load h5ad file
adata = sc.read_h5ad('sample_processed.h5ad')

# Continue analysis
sc.tl.rank_genes_groups(adata, 'seurat_clusters', method='wilcoxon')
```

#### Continue in R/Seurat

```r
library(Seurat)

# Load processed object
seurat_obj <- readRDS('sample_processed.rds')

# Continue analysis
seurat_obj <- FindSubCluster(seurat_obj, cluster = 0, resolution = 0.5)
```

## Getting Help

- **Documentation**: Check docs/ directory
- **Issues**: https://github.com/salzcamino/nf-scrnaseq/issues
- **Nextflow help**: https://www.nextflow.io/docs/latest/
