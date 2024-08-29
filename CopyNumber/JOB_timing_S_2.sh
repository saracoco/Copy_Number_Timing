#!/bin/bash
#SBATCH --job-name=Timing
#SBATCH --no-requeue
#SBATCH --nodes=2
#SBATCH --cpus-per-task=24
#SBATCH --time=10:00:00
#SBATCH --partition=THIN
#SBATCH --mem=400gb
#SBATCH --output=SN_3.out

module load R

R CMD BATCH scripts/simulate_S_2.R


module purge
