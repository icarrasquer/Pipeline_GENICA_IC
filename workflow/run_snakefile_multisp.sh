#!/bin/bash
#SBATCH --job-name=cleaning_mapping_pipeline
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=400G
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu_long
#SBATCH --account=paygo
#SBATCH --wckey=ips_genica
#SBATCH --output=jobs/smk_%j.out
#SBATCH --error=jobs/smk_%j.err

module load Anaconda3
source $(conda info --base)/etc/profile.d/conda.sh

conda activate snakemake

snakemake --unlock --snakefile Snakefile_multisp.smk
  
snakemake \
  --snakefile Snakefile_multisp.smk \
  --configfile config/config_Piun.yaml \
  --use-conda --rerun-incomplete  --rerun-triggers mtime \
  --cores ${SLURM_CPUS_PER_TASK} \
  --resources repair_slots=1 markdup_slots=2 mem_mb=70000
