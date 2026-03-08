#!/bin/bash
#SBATCH --job-name=read_cleaning
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --account=paygo
#SBATCH --wckey=ips_genica
#SBATCH --output=jobs/smk_%j.out
#SBATCH --error=jobs/smk_%j.err

module load Anaconda3
source $(conda info --base)/etc/profile.d/conda.sh
# activate the env that has snakemake (and maybe mamba)
conda activate snakemake

snakemake --unlock --snakefile workflow/Snakefile_droc_oldlibprep.smk

snakemake \
  --snakefile workflow/Snakefile_droc_oldlibprep.smk \
  --configfile config/config_droc_oldlibprep.yaml \
  --use-conda --rerun-incomplete  --rerun-triggers mtime \
  --cores ${SLURM_CPUS_PER_TASK} \
  --resources repair_slots=1 mem_mb=64000
