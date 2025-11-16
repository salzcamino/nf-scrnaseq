# Claude Code Guide for nf-scrnaseq

## Project Overview
Nextflow DSL2 pipeline for single-cell RNA sequencing analysis using Scanpy. Processes scRNA-seq data through QC, integration, clustering, annotation, and functional analysis.

## Architecture

```
main.nf                     # Main workflow orchestrator
modules/local/*.nf          # 17 process modules (Python embedded in Nextflow)
nextflow.config             # Global config with 60+ parameters
conf/*.config               # Cloud/HPC execution profiles (AWS, GCP, SLURM, PBS)
environment.yml             # Conda dependencies
Dockerfile                  # Container definition
```

### Pipeline Flow
Import → QC → Doublet Detection → Integration → Normalize → Cell Cycle → HVG → Dim Reduction → Batch Correction → Clustering → Diff Expression → Annotation → GSEA → Trajectory → Cell Communication → Report → Export

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

## Common Tasks

### Add New Parameter
1. Add to `params` block in [nextflow.config](nextflow.config)
2. Add validation in main.nf if needed
3. Use in relevant module
4. Document in README.md parameter table

### Add New Integration Method
Modify [sample_integration.nf](modules/local/sample_integration.nf):
```python
elif method == 'new_method':
    # Implementation
```

### Add New Clustering Algorithm
Modify [clustering.nf](modules/local/clustering.nf):
```python
if clustering_algorithm == 'new_algo':
    # Implementation
```

### Add New Export Format
Add to [data_export.nf](modules/local/data_export.nf)

## Important Conventions

### Code Style
- Python code embedded in triple-quoted strings within Nextflow
- Use raw strings (no indentation issues) - fixed in recent commits
- Scanpy for all single-cell operations
- AnnData (.h5ad) as primary data format

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

## Quick Reference

### Run Full Pipeline
```bash
nextflow run main.nf --input_dir ./data --outdir ./results -profile conda
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
