#!/usr/bin/env python3
"""
Generate synthetic 10X Genomics format test data for nf-scrnaseq pipeline
"""

import random
from pathlib import Path

# Set random seed for reproducibility
random.seed(42)

# Parameters
n_cells = 50
n_genes = 100
output_dir = Path("10x_sample")

# Create output directory
output_dir.mkdir(exist_ok=True)

# Generate gene names (including MT and ribosomal genes)
gene_names = []
gene_ids = []

# Add some mitochondrial genes (10%)
n_mt_genes = 10
for i in range(n_mt_genes):
    gene_id = f"ENSG0000000{i:03d}"
    gene_name = f"MT-CO{i+1}" if i < 5 else f"MT-ND{i-4}"
    gene_names.append(gene_name)
    gene_ids.append(gene_id)

# Add some ribosomal genes (10%)
n_ribo_genes = 10
for i in range(n_ribo_genes):
    gene_id = f"ENSG0000001{i:03d}"
    gene_name = f"RPS{i+1}" if i < 5 else f"RPL{i-4}"
    gene_names.append(gene_name)
    gene_ids.append(gene_id)

# Add regular genes (80%)
n_regular_genes = n_genes - n_mt_genes - n_ribo_genes
for i in range(n_regular_genes):
    gene_id = f"ENSG0000002{i:03d}"
    gene_name = f"GENE{i+1}"
    gene_names.append(gene_name)
    gene_ids.append(gene_id)

# Generate cell barcodes
barcodes = [f"CELL{i+1:04d}-1" for i in range(n_cells)]

# Generate sparse count matrix
# Most entries will be zero (sparse), but we'll create some realistic patterns
matrix_data = []

for cell_idx in range(n_cells):
    # Each cell expresses a subset of genes
    n_expressed = random.randint(20, 60)  # Each cell expresses 20-60 genes
    expressed_genes = random.sample(range(n_genes), n_expressed)

    for gene_idx in expressed_genes:
        # Generate counts (higher for regular genes, lower for MT/ribo)
        if gene_idx < n_mt_genes:
            # MT genes: lower counts but present
            count = max(1, int(random.expovariate(1/5)))
        elif gene_idx < n_mt_genes + n_ribo_genes:
            # Ribosomal genes: moderate counts
            count = max(1, int(random.expovariate(1/10)))
        else:
            # Regular genes: variable counts
            count = max(1, int(random.expovariate(1/20)))

        if count > 0:
            # Store as (gene_idx, cell_idx, count) - 1-indexed for MTX format
            matrix_data.append((gene_idx + 1, cell_idx + 1, count))

# Write matrix.mtx file in Matrix Market format
with open(output_dir / "matrix.mtx", "w") as f:
    f.write("%%MatrixMarket matrix coordinate integer general\n")
    f.write("%\n")
    f.write(f"{n_genes} {n_cells} {len(matrix_data)}\n")

    for gene_idx, cell_idx, count in matrix_data:
        f.write(f"{gene_idx} {cell_idx} {count}\n")

# Write genes.tsv file
with open(output_dir / "genes.tsv", "w") as f:
    for gene_id, gene_name in zip(gene_ids, gene_names):
        f.write(f"{gene_id}\t{gene_name}\n")

# Write barcodes.tsv file
with open(output_dir / "barcodes.tsv", "w") as f:
    for barcode in barcodes:
        f.write(f"{barcode}\n")

print(f"Test data generated successfully!")
print(f"Location: {output_dir}")
print(f"Cells: {n_cells}")
print(f"Genes: {n_genes}")
print(f"  - Mitochondrial: {n_mt_genes}")
print(f"  - Ribosomal: {n_ribo_genes}")
print(f"  - Regular: {n_regular_genes}")
print(f"Non-zero entries: {len(matrix_data)}")
