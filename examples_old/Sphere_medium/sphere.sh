#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --partition=RT
#SBATCH --job-name=sphere
#SBATCH --comment="-"

# export OMP_NUM_THREADS=1
srun ../../src/cli.o sphere.cfg
