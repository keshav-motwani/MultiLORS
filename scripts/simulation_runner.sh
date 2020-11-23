#!/bin/bash

#SBATCH --array=1-16500
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=LORS
#SBATCH --mem-per-cpu=6gb
#SBATCH -t 1440

export PATH=${HOME}/miniconda3/envs/r_env/bin:$PATH

Rscript scripts/simulation_runner.R
