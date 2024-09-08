#!/bin/bash
#SBATCH --job-name=Timing
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=10:00:00
#SBATCH --partition=THIN
#SBATCH --mem=400gb
#SBATCH --output=simulate.out


module load R/4.3.3

R CMD BATCH scripts/simulate.R

module purge
