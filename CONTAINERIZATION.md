# Containerization Guide - nf-scrnaseq Pipeline

This document provides quick reference for containerization setup and usage.

## Quick Start

### Option 1: Docker (Recommended for Local Development)

```bash
# Build the image (one-time setup)
docker build -t nf-scrnaseq:latest .

# Run pipeline with Docker
nextflow run main.nf --input data/ -profile docker
```

### Option 2: Singularity (Recommended for HPC)

```bash
# Option A: Use Galaxy pre-built image (easiest)
nextflow run main.nf --input data/ -profile singularity

# Option B: Build from Docker image
singularity build nf-scrnaseq.sif docker://nf-scrnaseq:latest
nextflow run main.nf --input data/ -profile singularity_local
```

### Option 3: Conda (No Containers Required)

```bash
# Conda environment created automatically
nextflow run main.nf --input data/ -profile conda
```

## System Requirements

| Environment | Requirements |
|-------------|--------------|
| Docker | Docker >= 20.10, 2GB disk space |
| Singularity | Singularity >= 3.0, 2GB disk space |
| Conda | Conda or Miniconda, 3-4GB disk space |

## Complete Dependency List

The Dockerfile includes all Python packages from `environment.yml`:

### Core Analysis
- scanpy >= 1.9.3 - Single-cell analysis framework
- anndata >= 0.9.2 - Data structure for single-cell
- leidenalg >= 0.10.0 - Leiden clustering algorithm

### Data Science
- pandas >= 2.0 - Data manipulation
- numpy >= 1.24 - Numerical computing
- scipy >= 1.11 - Scientific computing
- scikit-learn >= 1.3 - Machine learning

### Visualization
- matplotlib-base >= 3.7 - Plotting library
- seaborn >= 0.12 - Statistical visualization

### Integration Methods
- harmonypy >= 0.0.9 - Harmony integration
- scanorama >= 1.7.0 - Scanorama integration
- bbknn >= 1.6.0 - Batch-balanced k-NN
- *(Optional) scvi-tools - scVI integration*

### Specialized Tools
- scrublet >= 0.2.3 - Doublet detection
- celltypist >= 1.6.0 - Cell type annotation
- gseapy >= 1.0.0 - Gene set enrichment
- cellphonedb >= 5.0.0 - Cell-cell communication

### File I/O & Dependencies
- h5py >= 3.8 - HDF5 file support
- pytables >= 3.7 - Table support
- pyyaml >= 6.0 - YAML parsing
- python-igraph - Graph library for Leiden

## Building Custom Images

### Modify Dependencies

1. Edit `Dockerfile` and add/remove packages in the appropriate RUN section
2. Rebuild:
```bash
docker build -t nf-scrnaseq:custom .
```

3. Update `nextflow.config`:
```groovy
docker {
    process.container = 'nf-scrnaseq:custom'
}
```

### Add System Packages

Edit Dockerfile's apt-get section:
```dockerfile
RUN apt-get install -y \
    new-system-package \
    another-package
```

## Profile Configuration

### Docker Profile Settings
```groovy
docker {
    docker.enabled = true
    docker.runOptions = '-u $(id -u):$(id -g)'  # User ID mapping
    docker.temp = '/tmp'                         # Temp directory
    docker.mountFlags = 'z'                      # SELinux support
    process.container = 'nf-scrnaseq:latest'
}
```

### Singularity Profile Settings
```groovy
singularity {
    singularity.enabled = true
    singularity.autoMounts = true                # Auto-mount directories
    singularity.cacheDir = '~/.singularity-cache'
    singularity.pullTimeout = '20 min'
    process.container = 'https://depot.galaxyproject.org/singularity/scanpy:1.9.3--pyhdfd78af_0'
}
```

### Conda Profile Settings
```groovy
conda {
    conda.enabled = true
    conda.useMamba = false
    conda.cacheDir = '~/.nextflow-conda-cache'
    process.conda = '${projectDir}/environment.yml'
}
```

## Common Commands

### Docker

```bash
# Build image
docker build -t nf-scrnaseq:latest .

# Run with test data
nextflow run main.nf -profile test,docker

# Run with custom data
nextflow run main.nf --input /data/samples -profile docker

# Push to registry
docker login
docker tag nf-scrnaseq:latest username/nf-scrnaseq:latest
docker push username/nf-scrnaseq:latest

# Clean up unused images
docker system prune -a

# View image size
docker images nf-scrnaseq
```

### Singularity

```bash
# Build from Docker image
singularity build nf-scrnaseq.sif docker://nf-scrnaseq:latest

# Run with Galaxy image
nextflow run main.nf --input /data/samples -profile singularity

# Run with local image
nextflow run main.nf --input /data/samples -profile singularity_local

# View image info
singularity inspect nf-scrnaseq.sif

# Clean cache
rm -rf ~/.singularity-cache
```

### Conda

```bash
# Run with Conda (builds environment first run)
nextflow run main.nf --input /data/samples -profile conda

# Set cache directory
export NXF_CONDA_CACHEDIR=/custom/path
nextflow run main.nf -profile conda

# Clean cache
rm -rf ~/.nextflow-conda-cache
```

## Troubleshooting

### Docker Issues

**Cannot connect to daemon**
```bash
# Linux: Start Docker daemon
sudo systemctl start docker

# macOS: Open Docker app
open -a Docker
```

**Permission denied**
```bash
# Option 1: Add user to group
sudo usermod -aG docker $USER
newgrp docker

# Option 2: Use sudo
sudo nextflow run main.nf -profile docker
```

**Insufficient disk space**
```bash
docker system prune -a  # Remove all unused images
```

### Singularity Issues

**Cannot build from Docker**
```bash
# Requires Docker daemon running
docker run hello-world  # Test Docker access

# Or use Apptainer (newer versions)
apptainer build nf-scrnaseq.sif docker://nf-scrnaseq:latest
```

**Timeout pulling images**
```bash
# Increase timeout in nextflow.config
singularity.pullTimeout = '30 min'
```

### Conda Issues

**Environment creation fails**
```bash
# Clear cache and retry
rm -rf ~/.nextflow-conda-cache
nextflow run main.nf -profile conda

# Or set custom location
export NXF_CONDA_CACHEDIR=/path/with/space
```

**Very slow first run**
```bash
# Install mamba for faster dependency resolution
conda install -c conda-forge mamba

# Enable in nextflow.config
conda.useMamba = true  # Or set to false if it causes issues
```

## Performance Tips

1. **Docker**: Use built image cache - rebuilds are fast if no changes
2. **Singularity**: Pull images once, cache directory stores them
3. **Conda**: First run is slow; subsequent runs reuse environment
4. **General**: Use `-qs 20` to limit Nextflow queue size and avoid system overload

## Registry Options

### Docker Hub
```bash
docker push username/nf-scrnaseq:latest
```
Update nextflow.config:
```groovy
process.container = 'username/nf-scrnaseq:latest'
```

### Quay.io
```bash
docker push quay.io/yourorg/nf-scrnaseq:latest
```
Update nextflow.config:
```groovy
docker.registry = 'quay.io'
process.container = 'quay.io/yourorg/nf-scrnaseq:latest'
```

### Private Registry
```bash
docker push private.registry.com/nf-scrnaseq:latest
```
Update nextflow.config:
```groovy
docker.registry = 'private.registry.com'
process.container = 'private.registry.com/nf-scrnaseq:latest'
```

## CI/CD Integration

### GitHub Actions Example
```yaml
- name: Build Docker Image
  run: docker build -t nf-scrnaseq:test .

- name: Run Tests
  run: nextflow run main.nf -profile test,docker
```

### GitLab CI Example
```yaml
test:docker:
  image: docker:latest
  script:
    - docker build -t nf-scrnaseq:test .
    - docker run nf-scrnaseq:test nextflow run main.nf -profile test,docker
```

## See Also

- Full documentation: `README.md` - Containerization section
- Config reference: `nextflow.config` - Profiles section
- Package list: `environment.yml` - All conda dependencies
- Dockerfile: `Dockerfile` - Container definition

