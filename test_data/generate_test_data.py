#!/usr/bin/env python3
"""
Generate realistic 10X Genomics format test data for nf-scrnaseq pipeline.

Uses real gene symbols from human PBMC cells to enable proper testing of:
- Cell type annotation (CellTypist, marker scoring)
- Differential expression analysis
- Clustering and visualization

Creates synthetic expression data with realistic patterns for:
- T cells, B cells, Monocytes, NK cells
"""

import random
import math
from pathlib import Path
from collections import Counter

# Set random seed for reproducibility
random.seed(42)

# Parameters
n_cells = 200  # More cells for better clustering
output_dir = Path("10x_sample")

# Create output directory
output_dir.mkdir(exist_ok=True)

# Define real genes by category
# These are actual human gene symbols recognized by CellTypist and standard annotation tools

# Mitochondrial genes (for QC)
mt_genes = [
    "MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2",
    "MT-ND3", "MT-ND4", "MT-ND5", "MT-ATP6", "MT-CYB"
]

# Ribosomal genes (for QC)
ribo_genes = [
    "RPS2", "RPS3", "RPS4X", "RPS5", "RPS6",
    "RPL3", "RPL4", "RPL5", "RPL6", "RPL7"
]

# T cell markers
t_cell_genes = [
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B",
    "TRAC", "TRBC1", "TRBC2", "IL7R", "CCR7", "LEF1"
]

# B cell markers
b_cell_genes = [
    "CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "BANK1",
    "BLK", "CD22", "CD24", "IGHM", "IGHD", "TCL1A"
]

# Monocyte/Macrophage markers
monocyte_genes = [
    "CD14", "LYZ", "FCGR3A", "MS4A7", "CSF1R", "CD68",
    "ITGAM", "S100A8", "S100A9", "VCAN", "FCN1", "CTSS"
]

# NK cell markers
nk_cell_genes = [
    "NKG7", "GNLY", "KLRD1", "KLRF1", "KLRB1", "NCAM1",
    "FCGR3A", "PRF1", "GZMA", "GZMB", "KLRC1", "KLRK1"
]

# Dendritic cell markers
dc_genes = [
    "FCER1A", "CLEC10A", "CD1C", "ITGAX", "HLA-DRA",
    "HLA-DRB1", "HLA-DQA1", "HLA-DQB1"
]

# Platelet markers
platelet_genes = [
    "PPBP", "PF4", "GP9", "ITGA2B", "SELP", "CD9"
]

# Housekeeping/ubiquitous genes
housekeeping_genes = [
    "ACTB", "GAPDH", "B2M", "MALAT1", "TMSB4X", "EEF1A1",
    "TPT1", "UBC", "HSP90AA1", "HNRNPA1", "PTMA", "CALM1",
    "FTL", "FTH1", "S100A4", "LGALS1", "ANXA1", "ANXA2",
    "LDHA", "PGK1", "ENO1", "PKM", "VIM", "HSPA8"
]

# Combine all genes
all_genes = (
    mt_genes + ribo_genes + t_cell_genes + b_cell_genes +
    monocyte_genes + nk_cell_genes + dc_genes + platelet_genes +
    housekeeping_genes
)

# Remove duplicates while preserving order
seen = set()
gene_names = []
for gene in all_genes:
    if gene not in seen:
        gene_names.append(gene)
        seen.add(gene)

n_genes = len(gene_names)

# Generate synthetic Ensembl IDs
gene_ids = [f"ENSG{i:011d}" for i in range(n_genes)]

# Generate cell barcodes
barcodes = [f"CELL{i+1:04d}-1" for i in range(n_cells)]

# Define cell type composition (similar to PBMC)
# Assign each cell to a type with realistic proportions
cell_types = []
type_proportions = {
    'T_cell': 0.50,      # 50% T cells (most common in PBMC)
    'B_cell': 0.10,      # 10% B cells
    'Monocyte': 0.20,    # 20% Monocytes
    'NK_cell': 0.15,     # 15% NK cells
    'DC': 0.03,          # 3% Dendritic cells
    'Platelet': 0.02     # 2% Platelets
}

for _ in range(n_cells):
    rand_val = random.random()
    cumsum = 0
    for cell_type, prop in type_proportions.items():
        cumsum += prop
        if rand_val <= cumsum:
            cell_types.append(cell_type)
            break

# Define which genes are markers for each cell type
cell_type_markers = {
    'T_cell': set(t_cell_genes),
    'B_cell': set(b_cell_genes),
    'Monocyte': set(monocyte_genes),
    'NK_cell': set(nk_cell_genes),
    'DC': set(dc_genes),
    'Platelet': set(platelet_genes)
}

# Helper function for exponential distribution (without numpy)
def random_exponential(scale):
    """Generate exponential random variable with given scale (mean)."""
    return -scale * math.log(1.0 - random.random())

# Generate sparse count matrix
matrix_data = []

print("Generating expression matrix...")
for cell_idx in range(n_cells):
    cell_type = cell_types[cell_idx]

    for gene_idx, gene_name in enumerate(gene_names):
        # Base expression probability
        if gene_name in mt_genes:
            # MT genes: expressed in all cells at low-moderate levels
            if random.random() < 0.8:
                count = max(1, int(random_exponential(15)))
            else:
                continue
        elif gene_name in ribo_genes:
            # Ribosomal genes: highly expressed in all cells
            if random.random() < 0.9:
                count = max(1, int(random_exponential(50)))
            else:
                continue
        elif gene_name in housekeeping_genes:
            # Housekeeping: expressed in all cells
            if random.random() < 0.85:
                count = max(1, int(random_exponential(30)))
            else:
                continue
        elif gene_name in cell_type_markers.get(cell_type, set()):
            # Marker gene for this cell type: HIGH expression
            if random.random() < 0.85:
                count = max(1, int(random_exponential(40)))
            else:
                continue
        else:
            # Non-marker gene: LOW or no expression
            if random.random() < 0.15:  # 15% chance of background expression
                count = max(1, int(random_exponential(5)))
            else:
                continue

        # Add noise and dropout
        if count > 0 and random.random() > 0.1:  # 10% dropout
            # Store as (gene_idx, cell_idx, count) - 1-indexed for MTX format
            matrix_data.append((gene_idx + 1, cell_idx + 1, count))

print(f"Generated {len(matrix_data)} non-zero entries")

# Write matrix.mtx file in Matrix Market format
print("Writing matrix.mtx...")
with open(output_dir / "matrix.mtx", "w") as f:
    f.write("%%MatrixMarket matrix coordinate integer general\n")
    f.write("%\n")
    f.write(f"{n_genes} {n_cells} {len(matrix_data)}\n")

    for gene_idx, cell_idx, count in matrix_data:
        f.write(f"{gene_idx} {cell_idx} {count}\n")

# Write genes.tsv file (standard 10X format)
print("Writing genes.tsv...")
with open(output_dir / "genes.tsv", "w") as f:
    for gene_id, gene_name in zip(gene_ids, gene_names):
        f.write(f"{gene_id}\t{gene_name}\n")

# Write barcodes.tsv file
print("Writing barcodes.tsv...")
with open(output_dir / "barcodes.tsv", "w") as f:
    for barcode in barcodes:
        f.write(f"{barcode}\n")

# Write cell type ground truth (for validation)
print("Writing cell_types_ground_truth.csv...")
with open(output_dir / "cell_types_ground_truth.csv", "w") as f:
    f.write("barcode,cell_type\n")
    for barcode, cell_type in zip(barcodes, cell_types):
        f.write(f"{barcode},{cell_type}\n")

# Summary statistics
print(f"\n" + "="*50)
print(f"Test data generated successfully!")
print(f"="*50)
print(f"Location: {output_dir}")
print(f"Cells: {n_cells}")
print(f"Genes: {n_genes}")
print(f"  - Mitochondrial: {len(mt_genes)}")
print(f"  - Ribosomal: {len(ribo_genes)}")
print(f"  - T cell markers: {len(t_cell_genes)}")
print(f"  - B cell markers: {len(b_cell_genes)}")
print(f"  - Monocyte markers: {len(monocyte_genes)}")
print(f"  - NK cell markers: {len(nk_cell_genes)}")
print(f"  - DC markers: {len(dc_genes)}")
print(f"  - Platelet markers: {len(platelet_genes)}")
print(f"  - Housekeeping: {len(housekeeping_genes)}")
print(f"Non-zero entries: {len(matrix_data)}")
print(f"\nCell type composition:")
type_counts = Counter(cell_types)
for cell_type, count in sorted(type_counts.items()):
    print(f"  {cell_type}: {count} cells ({100*count/n_cells:.1f}%)")
print(f"\nGround truth saved to: {output_dir / 'cell_types_ground_truth.csv'}")
