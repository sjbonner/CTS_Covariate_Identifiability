#!/bin/bash
#SBATCH --array=1-240
#SBATCH --time=2:00:00
#SBATCH --mem=5G

module load r

Rscript --vanilla run_replicate.R $SLURM_ARRAY_TASK_ID &> Logs/rep_${SLURM_ARRAY_TASK_ID}.Rout
