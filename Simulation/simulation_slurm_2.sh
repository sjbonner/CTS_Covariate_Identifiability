#!/bin/bash
#SBATCH --array=21-50,71-100,121-150,171-200,221-250,271-300,321-350,371-400,421-450
#SBATCH --time=4:00:00
#SBATCH --mem=5G
#SBATCH --error=Logs/rep_${SLURM_ARRAY_TASK_ID}.Rout
#SBATCH --output=Logs/rep_${SLURM_ARRAY_TASK_ID}.Rout


module load gcc/7.3.0
module load r

Rscript --vanilla run_replicate.R $SLURM_ARRAY_TASK_ID &> Logs/rep_${SLURM_ARRAY_TASK_ID}.Rout
