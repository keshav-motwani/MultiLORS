#!/bin/bash

#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=MultiLORS_tune
#SBATCH --mem-per-cpu=6gb
#SBATCH -t 4320
#SBATCH --output=logs/MultiLORS_tune.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --constraint=skylake

export OMP_NUM_THREADS=1

ml intel/2019
ml R/3.6.2

Rscript scripts/applications/hao_3_prime/MultiLORS_tune.R
