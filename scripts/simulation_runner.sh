#!/bin/bash

#SBATCH --array=1-3000%500
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=MultiLORS
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 1440
#SBATCH --output=logs/MultiLORS_%a.log

export PATH=${HOME}/miniconda3/envs/r_env/bin:$PATH

Rscript scripts/simulation_runner.R
