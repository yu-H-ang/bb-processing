#!/bin/sh
#
#SBATCH --partition=shared
#SBATCH --time=100:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=ZYH
#SBATCH --output=main.out
#SBATCH --open-mode=trnucate

srun main
