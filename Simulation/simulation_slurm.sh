#!/bin/bash
#SBATCH --array=1-5
#SBATCH --time=7:00:00
#SBATCH --mem=5G

module load gcc/7.3.0
module load r

Rscript --vanilla run_replicate.R $SLURM_ARRAY_TASK_ID &> Logs/rep_${SLURM_ARRAY_TASK_ID}.Rout
