process IMPORT_DATA {
    tag "Importing data"
    label 'process_low'

    publishDir "${params.outdir}/import", mode: params.publish_dir_mode

    input:
    path input_data
    val input_format

    output:
    path "raw_data.h5ad", emit: adata
    path "import_summary.txt", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3

    import sys
    import scanpy as sc
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import yaml

    # Set scanpy settings
    sc.settings.verbosity = 3

    # Determine input format
    input_path = Path('${input_data}')
    format_type = '${input_format}'

    print(f"Input path: {input_path}")
    print(f"Format type: {format_type}")
    print(f"Path exists: {input_path.exists()}")
    print(f"Is directory: {input_path.is_dir()}")

    # Auto-detect format if needed
    if format_type == 'auto':
        if input_path.is_dir():
            # Check if it's a 10X directory
            if (input_path / 'matrix.mtx.gz').exists() or (input_path / 'matrix.mtx').exists():
                format_type = '10x'
            else:
                raise ValueError(f"Cannot auto-detect format for directory: {input_path}")
        elif input_path.suffix == '.h5ad':
            format_type = 'h5ad'
        elif input_path.suffix in ['.csv', '.tsv', '.txt']:
            format_type = 'csv'
        elif input_path.suffix == '.h5':
            format_type = '10x_h5'
        else:
            raise ValueError(f"Cannot auto-detect format for file: {input_path}")

    print(f"Detected format: {format_type}")

    # Load data based on format
    if format_type == '10x':
        adata = sc.read_10x_mtx(
            input_path,
            var_names='gene_symbols',
            cache=False
        )
    elif format_type == '10x_h5':
        adata = sc.read_10x_h5(input_path)
    elif format_type == 'h5ad':
        adata = sc.read_h5ad(input_path)
    elif format_type == 'csv':
        # Assume genes as rows, cells as columns
        df = pd.read_csv(input_path, index_col=0)
        adata = sc.AnnData(df.T)
    else:
        raise ValueError(f"Unsupported format: {format_type}")

    # Basic statistics
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    total_counts = adata.X.sum()

    print(f"Loaded {n_cells} cells and {n_genes} genes")
    print(f"Total counts: {total_counts:,.0f}")

    # Save raw data
    adata.write('raw_data.h5ad', compression='gzip')

    # Write summary
    with open('import_summary.txt', 'w') as f:
        f.write(f"Data Import Summary\\n")
        f.write(f"==================\\n")
        f.write(f"Input format: {format_type}\\n")
        f.write(f"Number of cells: {n_cells:,}\\n")
        f.write(f"Number of genes: {n_genes:,}\\n")
        f.write(f"Total counts: {total_counts:,.0f}\\n")
        f.write(f"Mean counts per cell: {total_counts/n_cells:,.1f}\\n")
        f.write(f"Median counts per cell: {np.median(np.array(adata.X.sum(axis=1)).flatten()):,.1f}\\n")

    # Write versions
    versions = {
        'IMPORT_DATA': {
            'scanpy': sc.__version__,
            'pandas': pd.__version__,
            'numpy': np.__version__
        }
    }
    with open('versions.yml', 'w') as f:
        yaml.dump(versions, f)

    print("Data import completed successfully!")
    """

    stub:
    """
    touch raw_data.h5ad
    touch import_summary.txt
    echo "IMPORT_DATA:" > versions.yml
    """
}
