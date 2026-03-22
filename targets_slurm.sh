#!/bin/bash

#SBATCH --cpus-per-task 3
#SBATCH --mem 1000G

Rscript targets_slurm.R


#CONDA_ENV_NAME="leptin_paper"
#eval "$(conda shell.bash hook)"
#Activate the Conda environment
#conda activate "$CONDA_ENV_NAME"

# Check if Conda environment was successfully activated
#if [ $? -eq 0 ]; then
  #echo "Conda environment '$CONDA_ENV_NAME' activated."
#else
  #echo "Failed to activate Conda environment '$CONDA_ENV_NAME'."
#  exit 1
#fi
