# Changelog

All notable changes to the nf-scrnaseq project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive pytest test suite with 30+ tests
  - Module syntax validation
  - Workflow structure tests
  - Configuration validation
  - Test data integrity checks
  - Documentation completeness tests
- GitHub Actions CI/CD workflows
  - Automated testing on push and PR
  - Python 3.9, 3.10, 3.11 matrix testing
  - Code linting with flake8
  - Security scanning with bandit and safety
  - Nextflow syntax validation
- Release workflow for Docker image publishing
- Development dependencies file (`requirements-dev.txt`)
- Test suite documentation (`tests/README.md`)
- Functional test report (`FUNCTIONAL_TEST_REPORT.md`)
- This CHANGELOG file

### Fixed
- Verified exception handling in all modules (no issues found - multiple except clauses are intentional)

### Documentation
- Added comprehensive test coverage documentation
- Documented CI/CD pipeline setup
- Added instructions for running tests locally

## [0.1.0] - 2025-01-XX

### Added
- Initial release of nf-scrnaseq pipeline
- 17 analysis modules covering full scRNA-seq workflow
- Data import supporting 10X Genomics, H5AD, CSV formats
- Quality control and filtering
- Doublet detection (Scrublet, scDblFinder, DecontX)
- Multi-sample integration (Harmony, Scanorama, BBKNN, scVI)
- Normalization and log transformation
- Highly variable gene selection
- Dimensionality reduction (PCA, UMAP, t-SNE)
- Batch correction with automatic effect detection
- Cell cycle scoring and optional regression
- Clustering (Leiden, Louvain, Seurat SNN, Celda)
- Differential expression analysis
- Cell type annotation (CellTypist, marker-based)
- Gene set enrichment analysis (GSEA)
- Trajectory analysis (PAGA, diffusion pseudotime)
- Cell-cell communication analysis
- HTML report generation
- Multi-format data export (Seurat, Loom, CellxGene, CSV, MTX)
- Docker containerization
- Singularity support
- Conda environment specification
- Cloud execution profiles (AWS Batch, GCP, Slurm, PBS)
- Comprehensive documentation (README, guides)
- Realistic test dataset (200 cells, 105 genes with ground truth)

### Documentation
- Extensive README with parameter tables and examples
- Quick start guide
- Cloud execution guide
- Containerization guide
- Developer guide (CLAUDE.md)

## Release Notes

### Version 0.1.0

First production-ready release of nf-scrnaseq, a comprehensive Nextflow pipeline for single-cell RNA-seq analysis.

**Key Features:**
- Complete scRNA-seq workflow from import to export
- 17 modular analysis steps
- Multiple integration and batch correction methods
- Automated cell type annotation
- Cloud and HPC execution support
- Extensive documentation

**Tested On:**
- Python 3.9, 3.10, 3.11
- Nextflow 22.10.0+
- Docker, Singularity, Conda environments
- Local and cloud execution

**Known Limitations:**
- R-based tools (scDblFinder, DecontX, Seurat, Celda) require manual R package installation
- scVI integration method requires separate installation of scvi-tools
- Large datasets may require significant memory (32GB+ recommended)

---

## Upgrade Guide

### From Pre-release to 0.1.0

No breaking changes. This is the first official release.

### Adding CI/CD

If upgrading an existing installation:

1. Install development dependencies:
   ```bash
   pip install -r requirements-dev.txt
   ```

2. Run tests:
   ```bash
   pytest
   ```

3. Enable GitHub Actions:
   - Workflows in `.github/workflows/` will run automatically
   - No additional setup required

---

## Contributing

See testing guidelines in `tests/README.md` before submitting pull requests.

All notable changes should be documented in this file under "Unreleased" section.
