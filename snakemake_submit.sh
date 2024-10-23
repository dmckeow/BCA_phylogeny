#!/bin/bash
#
#SBATCH --job-name=snakemake_main_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=48:10:00
#SBATCH --mem-per-cpu=300M
#SBATCH --output=slurm_logs/%x-%j.log
#SBATCH --qos=pipeline

mkdir -p slurm_logs
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

date
srun snakemake --use-conda -j1 --profile=cubi-v1
date
