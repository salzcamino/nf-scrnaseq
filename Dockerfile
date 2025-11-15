FROM python:3.10-slim

LABEL maintainer="salzcamino@gmail.com"
LABEL description="Docker container for nf-scrnaseq pipeline with scanpy and doublet detection"

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        gcc \
        g++ \
        libhdf5-dev \
        pkg-config && \
    rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
    scanpy==1.9.6 \
    anndata==0.10.3 \
    pandas==2.1.3 \
    numpy==1.26.2 \
    scipy==1.11.4 \
    matplotlib==3.8.2 \
    seaborn==0.13.0 \
    h5py==3.10.0 \
    tables==3.9.2 \
    pyyaml==6.0.1 \
    leidenalg==0.10.2 \
    python-igraph==0.11.3 \
    scikit-learn==1.3.2 \
    umap-learn==0.5.5 \
    scrublet==0.2.3 \
    numba==0.58.1

# Note: R-based tools (scDblFinder, DecontX, SoupX) are available only in the conda profile
# This Docker image supports Python-based doublet detection (Scrublet) and ambient RNA estimation

# Set working directory
WORKDIR /work

# Default command
CMD ["/bin/bash"]
