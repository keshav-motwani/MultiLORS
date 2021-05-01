#!/bin/bash

#SBATCH --array=0-1649
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=MultiLORS
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1440
#SBATCH --output=logs/MultiLORS_%a.log

ml R

Rscript scripts/simulation_runner.R
