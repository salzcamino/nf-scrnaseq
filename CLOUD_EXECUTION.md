# Cloud Execution Setup Guide for nf-scrnaseq

This guide provides a quick reference for setting up and running the nf-scrnaseq pipeline on different cloud and HPC platforms.

## Files Created

### Configuration Files
1. **conf/aws.config** - AWS Batch execution profile
2. **conf/gcp.config** - Google Cloud Platform (Life Sciences and Batch) profiles
3. **conf/slurm.config** - Slurm HPC scheduler profile
4. **conf/pbs.config** - PBS/Torque HPC scheduler profile

All profiles are automatically included in the main `nextflow.config` via `includeConfig` directives.

## Quick Start

### AWS Batch

```bash
# Configure AWS
export AWS_REGION=us-east-1
aws configure

# Create compute environment and queue (one-time setup)
aws batch create-compute-environment \
  --compute-environment-name nf-scrnaseq-compute-env \
  --type MANAGED

# Run pipeline on AWS
nextflow run main.nf \
  --input s3://my-bucket/data.h5ad \
  -profile aws
```

**Key Parameters**:
- `aws.region`: AWS region (default: us-east-1)
- `aws.batch.computeEnvironment`: Compute environment name
- `aws.batch.jobQueue`: Job queue names for different resource levels
- `process.queue`: Queue assignment per process label

### Google Cloud Platform

```bash
# Authenticate with GCP
gcloud auth login
gcloud config set project MY_PROJECT_ID

# Create service account
gcloud iam service-accounts create nf-scrnaseq
gcloud projects add-iam-policy-binding MY_PROJECT_ID \
  --member=serviceAccount:nf-scrnaseq@MY_PROJECT_ID.iam.gserviceaccount.com \
  --role=roles/lifesciences.workflowsWorker

# Run pipeline on GCP
nextflow run main.nf \
  --input gs://my-bucket/data.h5ad \
  -profile gcp
```

**Key Parameters**:
- `google.project`: GCP project ID
- `google.location`: Compute region (e.g., us-central1)
- `google.lifeSciences.bootDiskSize`: Boot disk size (default: 100.GB)
- `google.lifeSciences.preemptible`: Use preemptible instances (default: false)

**Alternative**: Use `gcp_batch` profile for Google Cloud Batch instead of Life Sciences.

### Slurm HPC

```bash
# Load required modules
module load singularity nextflow

# Edit conf/slurm.config with your queue names
# Edit params.outdir to use /scratch or cluster storage

# Run pipeline on HPC
nextflow run main.nf \
  --input /path/to/data.h5ad \
  -profile slurm

# Or submit as batch job
sbatch << 'SLURM'
#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB

module load singularity nextflow
nextflow run main.nf --input data.h5ad -profile slurm
SLURM
```

**Key Parameters**:
- `process.queue`: Default queue name (check with `sinfo`)
- `executor.queueSize`: Number of jobs to queue
- `singularity.cacheDir`: Singularity cache directory
- `process.scratch`: Temporary scratch space location

### PBS/Torque HPC

```bash
# Load required modules
module load singularity nextflow

# Run pipeline on PBS
nextflow run main.nf \
  --input /path/to/data.h5ad \
  -profile pbs

# Or submit as batch job
qsub << 'PBS'
#!/bin/bash
#PBS -l walltime=168:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=16GB

module load singularity nextflow
nextflow run main.nf --input data.h5ad -profile pbs
PBS
```

**Key Parameters**:
- `process.queue`: Default queue name (check with `qstat`)
- `executor.queueSize`: Number of jobs to queue
- `singularity.cacheDir`: Singularity cache directory
- `process.scratch`: Temporary scratch space location

## Configuration Parameters Reference

### AWS Batch (conf/aws.config)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `aws.region` | `us-east-1` | AWS region |
| `aws.batch.computeEnvironment` | `nf-scrnaseq-compute-env` | Compute environment name |
| `aws.batch.jobQueue['default']` | `nf-scrnaseq-queue` | Queue for normal jobs |
| `aws.batch.jobQueue['high-memory']` | `nf-scrnaseq-highmem-queue` | Queue for memory-intensive jobs |
| `params.outdir` | `s3://my-bucket/results` | S3 bucket for results |
| Process Resources | See below | Per-label specifications |

**Process Resource Specifications**:
- `process_low`: 2 CPUs, 4 GB, 4h timeout
- `process_medium`: 4 CPUs, 16 GB, 12h timeout
- `process_high`: 8 CPUs, 32 GB, 24h timeout

### Google Cloud Platform (conf/gcp.config)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `google.project` | `my-gcp-project` | GCP project ID |
| `google.location` | `us-central1` | Compute region |
| `google.lifeSciences.bootDiskSize` | `100.GB` | Boot disk size |
| `google.lifeSciences.preemptible` | `false` | Use preemptible instances |
| `params.outdir` | `gs://my-bucket/results` | GCS bucket for results |

**Resource Specifications**:
- `process_low`: 2 CPUs, 4 GB, 4h, n1-standard-2
- `process_medium`: 4 CPUs, 16 GB, 12h, n1-standard-4
- `process_high`: 8 CPUs, 32 GB, 24h, n1-standard-8

### Slurm HPC (conf/slurm.config)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `process.executor` | `slurm` | Job scheduler |
| `process.queue['normal']` | `normal` | Default queue |
| `process.queue['long']` | `long` | Queue for long-running jobs |
| `executor.queueSize` | `20` | Max concurrent jobs |
| `executor.submitRateLimit` | `2 sec` | Rate limit for submissions |
| `executor.pollInterval` | `30 sec` | Job status check interval |
| `singularity.enabled` | `true` | Use Singularity containers |
| `process.scratch` | `/scratch/$USER` | Temporary storage |
| `params.outdir` | `/scratch/$USER/results` | Output directory |

**Queue Mapping**:
- `process_low`: `normal` queue
- `process_medium`: `normal` queue
- `process_high`: `long` queue

### PBS/Torque HPC (conf/pbs.config)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `process.executor` | `pbs` | Job scheduler |
| `process.queue['standard']` | `standard` | Default queue |
| `process.queue['long']` | `long` | Queue for long-running jobs |
| `executor.queueSize` | `20` | Max concurrent jobs |
| `executor.submitRateLimit` | `2 sec` | Rate limit for submissions |
| `executor.pollInterval` | `30 sec` | Job status check interval |
| `singularity.enabled` | `true` | Use Singularity containers |
| `process.scratch` | `$TMPDIR` | Temporary storage |
| `params.outdir` | `/scratch/$USER/results` | Output directory |

## Customization

### Changing Resource Allocations

Edit the relevant profile config file to adjust resources:

```groovy
process {
    withLabel: process_high {
        cpus   = { check_max( 16, 'cpus' ) }      // Increase CPUs
        memory = { check_max( 64.GB * task.attempt, 'memory' ) }  // Increase memory
        time   = { check_max( 48.h * task.attempt, 'time' ) }    // Increase time
    }
}
```

### Changing Queue Names

For HPC systems, update queue names to match your cluster:

```groovy
// In conf/slurm.config or conf/pbs.config
process {
    withLabel: process_medium {
        queue = 'your-queue-name'
    }
}
```

### Using Preemptible Instances (Cost Savings)

For GCP:
```groovy
// In conf/gcp.config
google {
    lifeSciences {
        preemptible = true  // 70% cost reduction (can be interrupted)
    }
}
```

For AWS Batch, enable spot instances when creating the compute environment:
```bash
aws batch create-compute-environment \
  --compute-resources spotIamFleetRole=arn:aws:iam::ACCOUNT:role/AmazonEC2SpotFleetRole
```

### Custom Container Images

Update the container image reference in your profile:

```groovy
// For AWS Batch
docker {
    enabled = true
    process.container = 'YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/nf-scrnaseq:latest'
}

// For GCP
docker {
    enabled = true
    process.container = 'us-central1-docker.pkg.dev/YOUR_PROJECT/nf-scrnaseq/nf-scrnaseq:latest'
}
```

## Cost Estimation

### AWS Batch

- **EC2 m5.xlarge** (4 CPU, 16 GB): ~$0.192/hour
- **S3 Storage**: $0.023/GB/month
- **Typical small run** (10 cells): ~2-4 hours = $0.38-0.77
- **Typical large run** (50K cells): ~24-48 hours = $4.6-9.2

### Google Cloud Platform

- **n1-standard-4**: ~$0.19/hour ($0.048 preemptible)
- **GCS Storage**: $0.020/GB/month
- **Typical small run**: ~2-4 hours = $0.38-0.76 ($0.10-0.19 preemptible)
- **Typical large run**: ~24-48 hours = $4.6-9.2 ($1.2-2.3 preemptible)

### HPC (Slurm/PBS)

- **No direct compute costs** (usually covered by institutional allocation)
- **Typical allocation**: 100,000 - 1,000,000 CPU hours per month
- **Storage quota**: Usually 1-10 TB

## Troubleshooting

### AWS Batch

**Problem**: Job failed to submit
- Check IAM permissions: `aws iam list-user-policies --user-name YOUR_USER`
- Verify compute environment: `aws batch describe-compute-environments`

**Problem**: Container image not found
- Verify ECR repository: `aws ecr describe-repositories`
- Check image tag: `aws ecr list-images --repository-name nf-scrnaseq`

### Google Cloud Platform

**Problem**: API not enabled
```bash
gcloud services enable lifesciences.googleapis.com
gcloud services enable compute.googleapis.com
gcloud services enable storage-api.googleapis.com
```

**Problem**: Authentication failed
```bash
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/key.json
gcloud auth application-default login
```

### HPC Clusters

**Problem**: Singularity image not found
- Check cache location: `ls -lh ~/.singularity-cache`
- Verify image permissions: `singularity exec /path/to/image.sif echo "OK"`

**Problem**: Job stays in queue
- Check queue status: `squeue` (Slurm) or `qstat` (PBS)
- Check resource availability: `sinfo` (Slurm) or `pbsnodes` (PBS)

## Additional Resources

- **Nextflow Documentation**: https://nextflow.io/docs
- **AWS Batch Guide**: https://docs.aws.amazon.com/batch/
- **GCP Life Sciences API**: https://cloud.google.com/life-sciences/docs
- **Slurm Documentation**: https://slurm.schedmd.com/
- **PBS Documentation**: https://www.pbspro.org/

## Support

For issues specific to cloud execution:
1. Check the main README.md Cloud Execution Profiles section
2. Review platform-specific documentation
3. Check Nextflow logs: `nextflow log`
4. Enable debug mode: `nextflow run main.nf -profile aws -debug`

