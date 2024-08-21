#!/bin/bash
#SBATCH --job-name=Timing
#SBATCH --no-requeue
#SBATCH --nodes=2
#SBATCH --cpus-per-task=24
#SBATCH --time=05:00:00
#SBATCH --partition=THIN
#SBATCH --mem=400gb

module load R

R CMD BATCH scripts/simulate.R


module purge
