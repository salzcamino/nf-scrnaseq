# Use Ubuntu base image for better system library support
FROM ubuntu:22.04

LABEL maintainer="salzcamino@gmail.com"
LABEL description="Docker container for nf-scrnaseq pipeline with scanpy, integration methods, and doublet detection"
LABEL version="1.0"

# Set environment variables to prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONNOUSERSITE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3.10 \
        python3-pip \
        python3-distutils \
        build-essential \
        gcc \
        g++ \
        gfortran \
        git \
        curl \
        wget \
        libhdf5-dev \
        libopenblas-dev \
        liblapack-dev \
        pkg-config \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

# Create Python symlink
RUN ln -s /usr/bin/python3.10 /usr/bin/python && \
    ln -s /usr/bin/pip3 /usr/bin/pip

# Upgrade pip, setuptools, and wheel
RUN pip install --upgrade \
    pip \
    setuptools \
    wheel

# Install core scientific Python packages
RUN pip install --no-cache-dir \
    numpy==1.26.2 \
    scipy==1.11.4 \
    pandas==2.1.3 \
    scikit-learn==1.3.2

# Install core single-cell analysis packages
RUN pip install --no-cache-dir \
    anndata==0.10.3 \
    scanpy==1.9.6 \
    leidenalg==0.10.2 \
    python-igraph==0.11.3 \
    umap-learn==0.5.5 \
    numba==0.58.1

# Install visualization packages
RUN pip install --no-cache-dir \
    matplotlib==3.8.2 \
    seaborn==0.13.0

# Install file I/O packages
RUN pip install --no-cache-dir \
    h5py==3.10.0 \
    tables==3.9.2 \
    pyyaml==6.0.1

# Install doublet detection (Python-based)
RUN pip install --no-cache-dir \
    scrublet==0.2.3

# Install integration methods
RUN pip install --no-cache-dir \
    harmonypy==0.0.9 \
    bbknn==1.6.0 \
    scanorama==1.7.0

# Install cell type annotation and communication tools
RUN pip install --no-cache-dir \
    celltypist==1.6.0 \
    cellphonedb==5.0.0

# Install gene set enrichment analysis
RUN pip install --no-cache-dir \
    gseapy==1.0.0

# Optional: scvi-tools for advanced integration
# Uncomment the following line if scVI integration is needed
# RUN pip install --no-cache-dir scvi-tools==1.0.1

# Note: R-based tools (scDblFinder, DecontX, SoupX) are NOT included in this Docker image
# These require R installation and bioconductor packages and are better handled via:
# 1. The conda profile (environment.yml)
# 2. A separate R container (multi-stage Docker build)
# Users requiring R-based tools should use the conda profile instead.

# Set working directory
WORKDIR /work

# Create directories for input/output
RUN mkdir -p /work/input /work/output

# Default command
CMD ["/bin/bash"]
