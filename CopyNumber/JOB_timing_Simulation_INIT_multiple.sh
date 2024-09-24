#!/bin/bash
#SBATCH --job-name=timing
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=15:00:00
#SBATCH --partition=EPYC
#SBATCH --mem=400gb
#SBATCH --output=simulate.out


module load R/4.3.3

R CMD BATCH scripts/simulate_with_initialization_multiple.R

module purge
