#!/bin/bash
#SBATCH --array=10-12
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH --error=Logs/rep_%A_%a_long.log
#SBATCH --output=Logs/rep_%A_%a_long.log

module load r

Rscript --vanilla run_replicate.R $SLURM_ARRAY_TASK_ID parameter_matrix_long.csv &> Logs/long_run_${SLURM_ARRAY_TASK_ID}.Rout
