#!/bin/bash

#SBATCH --array=1-1000
#SBATCH --nodes=1
#SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
#SBATCH --cpus-per-task=1

## turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wangrui8@uw.edu

ml fhR/4.0.4-foss-2020b

R CMD BATCH --no-save  cluster2.R Rout2/example_slurm${SLURM_ARRAY_TASK_ID}.Rout