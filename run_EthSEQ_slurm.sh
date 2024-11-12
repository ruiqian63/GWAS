#!/bin/bash

#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=48:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=EthSEQ_slurm_folate # Job name
#SBATCH --output=EthSEQ_slurm_folate.%j.out # Name of the output file
#SBATCH --error=EthSEQ_slurm_folate.%j.error
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of threads/tasks on one node
#SBATCH --mem=512G

export R_X=4.1
Rscript stratify_folate.R
