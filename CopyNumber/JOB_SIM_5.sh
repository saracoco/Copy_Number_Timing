#!/bin/bash
#SBATCH --job-name=timing
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=20:00:00
#SBATCH --partition=EPYC
#SBATCH --mem=400gb
#SBATCH --output=simple_init_.out


module load R/4.3.3

R CMD BATCH scripts/simulate_CNtime_5.R

module purge
