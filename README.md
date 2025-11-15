# nf-scrnaseq

A Nextflow pipeline for single-cell RNA-seq analysis, built following nf-core best practices.

## Pipeline Overview

This pipeline performs comprehensive single-cell RNA-seq analysis including:

- **Quality Control**: FastQC and MultiQC reporting
- **Alignment**: STAR-based alignment optimized for scRNA-seq
- **Quantification**: UMI-based gene counting using STAR Solo
- **Cell Filtering**: Quality-based cell filtering using Seurat
- **Normalization**: Log-normalization or SCTransform
- **Dimensionality Reduction**: PCA and UMAP
- **Clustering**: Graph-based clustering
- **Marker Detection**: Identification of cluster-specific genes
- **Visualization**: Comprehensive QC and analysis plots

## Quick Start

### Prerequisites

- Nextflow (≥22.10.0)
- Docker, Singularity, Podman, or Conda

### Installation

```bash
# Clone the repository
git clone https://github.com/salzcamino/nf-scrnaseq.git
cd nf-scrnaseq

# Run the test profile
nextflow run main.nf -profile test,docker
```

## Usage

### Basic Usage

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --genome GRCh38 \
  --protocol 10x_3prime \
  -profile docker
```

### Input Samplesheet Format

The pipeline requires a CSV samplesheet with the following columns:

| Column    | Description                          | Required |
|-----------|--------------------------------------|----------|
| sample    | Sample identifier                    | Yes      |
| fastq_1   | Path to read 1 FASTQ file           | Yes      |
| fastq_2   | Path to read 2 FASTQ file           | No       |
| protocol  | scRNA-seq protocol (overrides --protocol) | No |

Example samplesheet:

```csv
sample,fastq_1,fastq_2,protocol
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,10x_3prime
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,10x_3prime
```

### Reference Genome

You can specify a reference genome in two ways:

1. **Pre-built genome** (recommended for standard organisms):
```bash
--genome GRCh38  # or GRCm39 for mouse
```

2. **Custom genome**:
```bash
--fasta /path/to/genome.fa \
--gtf /path/to/genes.gtf
```

### Supported Protocols

- `10x_3prime` - 10x Genomics 3' gene expression (default)
- `10x_5prime` - 10x Genomics 5' gene expression
- `dropseq` - Drop-seq
- `smartseq2` - Smart-seq2 (full-length)
- `smartseq3` - Smart-seq3
- `indrop` - inDrop

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to input samplesheet (CSV) |
| `--outdir` | Output directory for results |
| `--genome` or `--fasta` & `--gtf` | Reference genome |

### Optional Parameters

#### General Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--protocol` | `10x_3prime` | scRNA-seq protocol |
| `--aligner` | `star` | Alignment tool (star, salmon) |
| `--skip_fastqc` | `false` | Skip FastQC step |
| `--skip_multiqc` | `false` | Skip MultiQC step |

#### Cell Filtering Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_genes` | `200` | Minimum genes per cell |
| `--min_cells` | `3` | Minimum cells per gene |
| `--max_mito_pct` | `5` | Maximum mitochondrial percentage |
| `--min_counts` | `1000` | Minimum UMI counts per cell |
| `--max_counts` | `null` | Maximum UMI counts per cell |
| `--max_genes` | `null` | Maximum genes per cell |

#### Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--n_hvgs` | `2000` | Number of highly variable genes |
| `--n_pcs` | `50` | Number of principal components |
| `--resolution` | `0.8` | Clustering resolution |
| `--normalization_method` | `lognorm` | Normalization method (lognorm, sctransform) |
| `--skip_clustering` | `false` | Skip clustering analysis |
| `--skip_annotation` | `false` | Skip cell type annotation |

#### Resource Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `240.h` | Maximum time per process |

## Output Directory Structure

```
results/
├── fastqc/                     # FastQC reports
├── star/                       # STAR alignment outputs
│   ├── *.bam                  # Aligned BAM files
│   └── log/                   # Alignment logs
├── qc/                        # Cell QC plots and filtered data
│   ├── *_qc_pre_filter.pdf   # Pre-filtering QC plots
│   ├── *_qc_post_filter.pdf  # Post-filtering QC plots
│   └── *_filtering_stats.csv # Filtering statistics
├── analysis/                  # Downstream analysis results
│   ├── *_processed.rds       # Processed Seurat objects
│   ├── *_markers.csv         # Cluster marker genes
│   ├── *_umap.pdf            # UMAP visualizations
│   ├── *_pca.pdf             # PCA plots
│   └── *_summary.txt         # Analysis summary
├── multiqc/                   # MultiQC report
│   └── multiqc_report.html   # Aggregated QC report
└── pipeline_info/             # Pipeline execution info
    ├── execution_timeline.html
    ├── execution_report.html
    └── pipeline_dag.html
```

## Key Output Files

### QC and Filtering

- **`*_qc_pre_filter.pdf`**: Violin plots and scatter plots showing cell QC metrics before filtering
- **`*_qc_post_filter.pdf`**: QC metrics after filtering
- **`*_filtering_stats.csv`**: Summary statistics of filtering process
- **`*_filtered.rds`**: Filtered Seurat object (R data format)

### Analysis Results

- **`*_processed.rds`**: Fully processed Seurat object with clustering
- **`*_markers.csv`**: All marker genes for each cluster
- **`*_top_markers.csv`**: Top 10 marker genes per cluster
- **`*_umap.pdf`**: UMAP plots colored by cluster and QC metrics
- **`*_pca.pdf`**: PCA plots and elbow plot
- **`*_hvg.pdf`**: Highly variable genes plot
- **`*.h5ad`**: AnnData format for Scanpy compatibility (if available)

## Profiles

The pipeline comes with several execution profiles:

- `docker` - Use Docker containers
- `singularity` - Use Singularity containers
- `conda` - Use Conda environments
- `test` - Run with minimal test data
- `test_full` - Run with full test dataset

Example:
```bash
nextflow run main.nf -profile docker,test
```

## Example Workflows

### Standard 10x Genomics Analysis

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --genome GRCh38 \
  --protocol 10x_3prime \
  --min_genes 200 \
  --max_mito_pct 10 \
  --resolution 0.5 \
  -profile docker
```

### High-Quality Cell Selection

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --genome GRCh38 \
  --min_genes 500 \
  --min_counts 2000 \
  --max_mito_pct 5 \
  --max_counts 50000 \
  -profile docker
```

### SCTransform Normalization

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --genome GRCh38 \
  --normalization_method sctransform \
  -profile docker
```

## Advanced Features

### Resuming Pipeline Execution

Nextflow supports automatic resumption of failed pipelines:

```bash
nextflow run main.nf -resume
```

### Custom Configuration

You can provide custom configuration files:

```bash
nextflow run main.nf -c custom.config
```

### Cloud Execution

The pipeline is compatible with cloud executors:

```bash
# AWS Batch
nextflow run main.nf -profile docker -bucket-dir s3://my-bucket/work

# Google Cloud
nextflow run main.nf -profile docker -bucket-dir gs://my-bucket/work
```

## Troubleshooting

### Common Issues

1. **Out of memory errors**: Increase `--max_memory` or reduce dataset size
2. **STAR alignment fails**: Ensure sufficient memory (≥32GB recommended)
3. **No cells pass filtering**: Adjust QC thresholds (`--min_genes`, `--max_mito_pct`)

### Getting Help

For issues and questions:
- Check the [documentation](docs/)
- Open an [issue](https://github.com/salzcamino/nf-scrnaseq/issues)

## Pipeline Architecture

The pipeline is organized following nf-core standards:

```
nf-scrnaseq/
├── main.nf                    # Main pipeline script
├── nextflow.config           # Main configuration
├── conf/                     # Configuration files
│   ├── base.config          # Resource configurations
│   ├── modules.config       # Module-specific settings
│   └── test.config          # Test profile
├── modules/                  # Process definitions
│   └── local/               # Local modules
├── workflows/               # Workflow definitions
├── bin/                     # Utility scripts
├── assets/                  # Pipeline assets
├── docs/                    # Documentation
└── tests/                   # Test data
```

## Credits

This pipeline was developed using the [Nextflow DSL2](https://www.nextflow.io/) workflow framework and follows [nf-core](https://nf-co.re/) best practices.

### Tools Used

- [Nextflow](https://www.nextflow.io/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [STAR](https://github.com/alexdobin/STAR)
- [Seurat](https://satijalab.org/seurat/)

## License

This project is licensed under the MIT License.

## Citation

If you use this pipeline, please cite:

```
nf-scrnaseq: A Nextflow pipeline for single-cell RNA-seq analysis
https://github.com/salzcamino/nf-scrnaseq
```
