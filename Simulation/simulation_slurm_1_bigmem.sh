#!/bin/bash
#SBATCH --array=451-600
#SBATCH --time=10:00:00
#SBATCH --mem=15G
#SBATCH --error=Logs/rep_%A_%a.log
#SBATCH --output=Logs/rep_%A_%a.log


module load gcc/7.3.0
module load r

Rscript --vanilla run_replicate.R $SLURM_ARRAY_TASK_ID parameter_matrix.csv &> Logs/rep_${SLURM_ARRAY_TASK_ID}.Rout
