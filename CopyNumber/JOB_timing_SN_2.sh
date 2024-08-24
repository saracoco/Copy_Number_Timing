#!/bin/bash
#SBATCH --job-name=Timing
#SBATCH --no-requeue
#SBATCH --nodes=2
#SBATCH --cpus-per-task=24
#SBATCH --time=10:00:00
#SBATCH --partition=THIN
#SBATCH --mem=400gb

module load R

R CMD BATCH scripts/simulate_night_2.R


module purge
