#!/bin/bash
#SBATCH --array=1-450
#SBATCH --time=4:00:00
#SBATCH --mem=5G
#SBATCH --error=Logs/rep_%A_%a.Rout
#SBATCH --output=Logs/rep_%A_%a.Rout


module load gcc/7.3.0
module load r

Rscript --vanilla run_replicate.R $SLURM_ARRAY_TASK_ID parameter_matrix.csv &> Logs/rep_${SLURM_ARRAY_TASK_ID}.Rout
