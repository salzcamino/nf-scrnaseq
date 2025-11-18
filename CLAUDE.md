# Claude Code Guide for nf-scrnaseq

## Project Overview
Nextflow DSL2 pipeline for single-cell RNA sequencing analysis using Scanpy. Processes scRNA-seq data through QC, integration, clustering, annotation, and functional analysis.

**Status:** Production-ready with comprehensive automated testing and CI/CD

## Architecture

```
main.nf                     # Main workflow orchestrator
modules/local/*.nf          # 17 process modules (Python embedded in Nextflow)
nextflow.config             # Global config with 60+ parameters
conf/*.config               # Cloud/HPC execution profiles (AWS, GCP, SLURM, PBS)
environment.yml             # Conda dependencies
Dockerfile                  # Container definition
tests/                      # Automated test suite (30+ tests)
.github/workflows/          # CI/CD pipelines
```

### Pipeline Flow
Import → QC → Doublet Detection → Integration → Normalize → Cell Cycle → HVG → Dim Reduction → Batch Correction → Clustering → Diff Expression → Annotation → GSEA → Trajectory → Cell Communication → Report → Export

## Quality Assurance

### Automated Testing (NEW)
- **30+ pytest tests** covering all modules
- **100% module coverage** - All 17 modules validated
- **Multi-version testing** - Python 3.9, 3.10, 3.11
- **CI/CD integration** - Tests run on every push/PR

### Test Suite
```
tests/
├── test_modules.py       # 30+ automated tests
├── conftest.py           # Pytest configuration
├── README.md             # Test documentation
└── __init__.py           # Test package
```

Run tests: `pytest -v`

## Key Files

### Core Pipeline
- [main.nf](main.nf) - Workflow definition and parameter validation
- [nextflow.config](nextflow.config) - All parameters, profiles, and resource configs
- [modules/local/](modules/local/) - Process modules with embedded Python scripts

### Important Modules
- [import_data.nf](modules/local/import_data.nf) - Data input (10X, H5AD, CSV)
- [sample_integration.nf](modules/local/sample_integration.nf) - Multi-sample merge (Harmony, Scanorama, BBKNN, scVI)
- [clustering.nf](modules/local/clustering.nf) - Leiden, Louvain, Seurat, Celda
- [cell_type_annotation.nf](modules/local/cell_type_annotation.nf) - CellTypist + marker scoring
- [batch_correction.nf](modules/local/batch_correction.nf) - Harmony, ComBat, BBKNN
- [data_export.nf](modules/local/data_export.nf) - Multi-format output

## Development Patterns

### Module Structure
Each module follows this pattern:
```nextflow
process MODULE_NAME {
    label 'process_medium'           // Resource requirements
    container 'image'                // Container spec
    conda 'environment.yml'          // Conda env

    input:
        path(input_h5ad)
        val(param)

    output:
        path("output.h5ad"), emit: adata
        path("*.pdf"), emit: plots
        path("*_summary.txt"), emit: summary

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    # Python code here
    """
}
```

### Parameter Access
Parameters defined in `nextflow.config` accessible via `params.param_name`:
```nextflow
params.min_genes = 200
params.leiden_resolution = 0.5
params.annotation_method = 'celltypist'
```

### Adding New Analysis
1. Create module in `modules/local/new_module.nf`
2. Import in main.nf: `include { NEW_MODULE } from './modules/local/new_module'`
3. Add to workflow chain
4. Define parameters in nextflow.config
5. Update documentation in README.md

## Testing

### Quick Test
```bash
nextflow run main.nf -profile test,conda
```
Uses test data: [test_data/10x_sample/](test_data/10x_sample/) (200 cells × 105 genes)

### Docker Test
```bash
docker build -t nf-scrnaseq:latest .
nextflow run main.nf -profile test,docker
```

### Resume Failed Run
```bash
nextflow run main.nf -resume
```

## Testing & CI/CD

### Running Tests Locally
```bash
# Install development dependencies
pip install -r requirements-dev.txt

# Run all tests
pytest -v

# Run specific test class
pytest tests/test_modules.py::TestModuleSyntax -v

# Run with coverage
pytest --cov=modules --cov-report=html
```

### CI/CD Pipeline
Automated workflows run on every push/PR:
- ✅ **Testing** - All 30+ tests across Python 3.9, 3.10, 3.11
- ✅ **Linting** - Code quality with flake8
- ✅ **Security** - Vulnerability scanning with bandit, safety
- ✅ **Nextflow** - Syntax validation

### Creating a Release
```bash
# Tag release
git tag -a v0.2.0 -m "Release v0.2.0"
git push origin v0.2.0

# GitHub Actions automatically:
# - Builds Docker image
# - Pushes to Docker Hub
# - Creates GitHub release
```

## Common Tasks

### Add New Parameter
1. Add to `params` block in [nextflow.config](nextflow.config)
2. Add validation in main.nf if needed
3. Use in relevant module
4. Document in README.md parameter table
5. **Add test** in `tests/test_modules.py`

### Add New Integration Method
1. Modify [sample_integration.nf](modules/local/sample_integration.nf):
```python
elif method == 'new_method':
    # Implementation
```
2. **Add test case** to verify new method

### Add New Clustering Algorithm
1. Modify [clustering.nf](modules/local/clustering.nf):
```python
if clustering_algorithm == 'new_algo':
    # Implementation
```
2. **Add test** for new algorithm

### Add New Export Format
1. Add to [data_export.nf](modules/local/data_export.nf)
2. **Add test** to verify export works

## Important Conventions

### Code Style
- Python code embedded in triple-quoted strings within Nextflow
- Use raw strings (no indentation issues) - fixed in recent commits
- Scanpy for all single-cell operations
- AnnData (.h5ad) as primary data format
- **Exception handling**: Multiple `except` clauses per `try` block are intentional
  ```python
  try:
      import optional_package
  except ImportError:
      # Handle missing dependency
  except Exception as e:
      # Handle other errors
  ```

### Code Quality
- **Automated testing** - All changes should include tests
- **CI/CD validation** - Tests run on every commit
- **Security scanning** - Automated vulnerability detection
- **Multi-version support** - Test on Python 3.9, 3.10, 3.11

### Resource Labels
- `process_low` - Light tasks (QC plots)
- `process_medium` - Standard analysis (clustering, DE)
- `process_high` - Memory intensive (integration, batch correction)

### Output Structure
```
results/
├── {module_name}/
│   ├── *.h5ad          # Primary data
│   ├── *.pdf           # Plots
│   ├── *.csv           # Tables
│   └── *_summary.txt   # Text summaries
└── report/
    └── scrnaseq_report.html
```

## Dependencies

### Core Stack
- Python 3.9-3.10
- Scanpy 1.9.6, AnnData 0.10.3
- Nextflow 22.10.0+

### Key Packages
- Integration: harmonypy, scanorama, bbknn
- Annotation: celltypist
- Communication: cellphonedb
- Enrichment: gseapy
- Visualization: matplotlib, seaborn

### R Dependencies (optional)
- scDblFinder, DecontX/Celda, Seurat (via rpy2)

## Git Ignored (Local Only)
- `work/` - Nextflow intermediate cache
- `results/` - Pipeline outputs
- `.nextflow.log*` - Execution logs
- `.nextflow/` - Nextflow history

## Caveats

1. **Embedded Python** - Python code lives inside Nextflow `.nf` files, not separate `.py` files
2. **Indentation** - Shell/Python blocks must have consistent indentation (recent bug fixes)
3. **Container size** - Docker image is large (2GB+) due to dependencies
4. **R tools** - Some features (scDblFinder, Seurat clustering) require R setup
5. **Memory** - Integration and batch correction can be memory-intensive
6. **File paths** - Nextflow handles paths; don't hardcode absolute paths in modules

## Recent Improvements (2025-11-18)

### Automated Testing Suite ✅
- **30+ pytest tests** added
- **100% module coverage** achieved
- **Multi-version testing** (Python 3.9, 3.10, 3.11)
- Test suite documentation added

### CI/CD Pipeline ✅
- GitHub Actions workflows configured
- Automated testing on every push/PR
- Security vulnerability scanning
- Automated Docker image publishing on release

### Code Quality ✅
- Comprehensive assessment completed (93.2% pass rate)
- Exception handling verified (all correct)
- Security audit passed (0 vulnerabilities)
- Documentation improved

### New Files
- `tests/` - Complete test suite
- `.github/workflows/` - CI/CD pipelines
- `requirements-dev.txt` - Development dependencies
- `pytest.ini` - Test configuration
- `CHANGELOG.md` - Version history
- `FUNCTIONAL_TEST_REPORT.md` - Assessment report
- `IMPROVEMENTS_SUMMARY.md` - Implementation summary

## Quick Reference

### Run Full Pipeline
```bash
nextflow run main.nf --input ./data --outdir ./results -profile conda
```

### Run Tests
```bash
# Install dev dependencies
pip install -r requirements-dev.txt

# Run all tests
pytest -v

# Run with coverage
pytest --cov=modules --cov-report=html
```

### Key Parameters to Check
- `params.input_type` - auto, 10x_mtx, h5ad, csv
- `params.sample_ids` - auto or comma-separated list
- `params.integration_method` - concat, harmony, scanorama, bbknn, scvi
- `params.clustering_algorithm` - leiden, louvain, seurat_snn, celda
- `params.annotation_method` - celltypist, marker_scoring, combined
- `params.batch_correction_method` - harmony, combat, bbknn, scanorama

### Execution Profiles
- `test` - Test parameters (reduced dataset)
- `conda` - Conda environment
- `docker` - Docker container
- `singularity` - Singularity container
- `aws`, `gcp`, `slurm`, `pbs` - Cloud/HPC profiles

### Development Workflow
1. Make changes to modules
2. Run tests: `pytest -v`
3. Check CI passes on GitHub
4. Create PR
5. Merge after approval
