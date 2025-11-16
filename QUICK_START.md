# Quick Start Guide - Containerization

This guide will help you get started with the nf-scrnaseq pipeline using containers.

## Choose Your Execution Method

### Option 1: Docker (Recommended for Local Development)

**Pros**: Self-contained, portable, good for development and sharing
**Cons**: Requires Docker daemon, small performance overhead

#### Step 1: Build the Image
```bash
cd /path/to/nf-scrnaseq
docker build -t nf-scrnaseq:latest .
```

This takes ~5-10 minutes on first run (downloads all packages).

#### Step 2: Verify the Build
```bash
docker images | grep nf-scrnaseq
# Output: nf-scrnaseq  latest  abc123...  2.5GB  ...
```

#### Step 3: Run the Pipeline
```bash
# Test run with included test data
nextflow run main.nf -profile test,docker

# Run with your data
nextflow run main.nf --input /path/to/data -profile docker
```

### Option 2: Singularity (Recommended for HPC)

**Pros**: HPC-friendly, no daemon, minimal overhead
**Cons**: Requires Singularity installation

#### Step 1a: Use Pre-built Galaxy Image (Easiest)
```bash
# No build needed - pulls from Galaxy automatically
nextflow run main.nf --input /path/to/data -profile singularity
```

#### Step 1b: Build Custom Image from Docker
```bash
# If you have Docker installed:
singularity build nf-scrnaseq.sif docker://nf-scrnaseq:latest

# Or from local Docker image:
singularity build nf-scrnaseq.sif docker-daemon://nf-scrnaseq:latest
```

#### Step 2: Run with Local Image
```bash
nextflow run main.nf --input /path/to/data -profile singularity_local
```

### Option 3: Conda (No Containers)

**Pros**: Simplest setup, no containers needed
**Cons**: First run is slow (~20-30 min)

#### Step 1: Run with Conda Profile
```bash
# Conda environment created automatically
nextflow run main.nf --input /path/to/data -profile conda

# Second run is much faster (environment cached)
nextflow run main.nf --input /path/to/data -profile conda
```

#### Clear Conda Cache (if needed)
```bash
rm -rf ~/.nextflow-conda-cache
```

## Common Workflow Examples

### Running Test Pipeline
```bash
# Docker
nextflow run main.nf -profile test,docker

# Singularity
nextflow run main.nf -profile test,singularity

# Conda
nextflow run main.nf -profile test,conda
```

### Running with Your Data
```bash
# Docker
nextflow run main.nf \
  --input /path/to/data \
  --outdir ./results \
  -profile docker

# Singularity
nextflow run main.nf \
  --input /path/to/data \
  --outdir ./results \
  -profile singularity

# Conda
nextflow run main.nf \
  --input /path/to/data \
  --outdir ./results \
  -profile conda
```

### Running with Custom Parameters
```bash
nextflow run main.nf \
  --input /path/to/data \
  --integration_method harmony \
  --n_top_genes 3000 \
  --leiden_resolution 0.8 \
  --max_genes 5000 \
  -profile docker
```

### Running with HPC Job Scheduler
```bash
# SLURM
nextflow run main.nf \
  --input /path/to/data \
  -profile singularity \
  -executor slurm

# SGE
nextflow run main.nf \
  --input /path/to/data \
  -profile singularity \
  -executor sge
```

## Troubleshooting

### Docker: "Cannot connect to Docker daemon"
```bash
# Start Docker service
sudo systemctl start docker     # Linux
open -a Docker                   # macOS
```

### Docker: "Permission denied"
```bash
# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker
# Log out and back in for changes to take effect
```

### Docker: "Insufficient space"
```bash
# Clean unused Docker resources
docker system prune -a
```

### Singularity: "Cannot build from Docker"
```bash
# Ensure Docker is running first
docker run hello-world

# Then build Singularity image
singularity build nf-scrnaseq.sif docker://nf-scrnaseq:latest
```

### Conda: "Very slow environment creation"
```bash
# Install mamba for faster resolution
conda install -c conda-forge mamba

# Enable in nextflow.config (optional)
# conda.useMamba = true
```

## Sharing Your Pipeline

### Share Docker Image

#### Push to Docker Hub
```bash
docker login
docker tag nf-scrnaseq:latest username/nf-scrnaseq:latest
docker push username/nf-scrnaseq:latest

# Share with collaborators
# They can run: nextflow run main.nf -profile docker
# (automatically uses username/nf-scrnaseq:latest)
```

#### Push to Quay.io
```bash
docker login quay.io
docker tag nf-scrnaseq:latest quay.io/yourorg/nf-scrnaseq:latest
docker push quay.io/yourorg/nf-scrnaseq:latest
```

## Performance Comparison

| Method | Setup Time | Runtime Overhead | First Run | Subsequent Runs |
|--------|-----------|-----------------|-----------|-----------------|
| Docker | 5-10 min | ~5% | Normal | Normal |
| Singularity (Galaxy) | 2-5 min | ~2% | ~1 min pull | Fast |
| Singularity (Local) | 5-10 min | ~2% | Fast | Fast |
| Conda | 0 min | ~0% | 20-30 min | 30 sec |

## System Requirements

| Tool | Requirement |
|------|-------------|
| Docker | >= 20.10, 2GB disk space |
| Singularity | >= 3.0, 2GB disk space |
| Conda | 3-4GB disk space |

## Next Steps

1. Choose a containerization method above
2. Run the test dataset to verify setup
3. Customize parameters as needed for your data
4. Check results in `./results` directory
5. See `README.md` for full documentation

## Getting Help

- Docker: See CONTAINERIZATION.md "Docker Issues" section
- Singularity: See CONTAINERIZATION.md "Singularity Issues" section
- Conda: See CONTAINERIZATION.md "Conda Issues" section
- Pipeline: See README.md for full parameter documentation

## For Developers

To modify the container setup:

1. **Add Python packages**: Edit `Dockerfile`, rebuild: `docker build -t nf-scrnaseq:latest .`
2. **Add system packages**: Edit Dockerfile's `apt-get install` section
3. **Use custom image**: Update `nextflow.config` process.container setting
4. **For R tools**: Use conda profile instead (R packages in environment.yml)

See CONTAINERIZATION.md for detailed customization instructions.
