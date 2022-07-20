#!/bin/bash

#SBATCH --job-name=aem1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=5-00:00:00
#SBATCH --mem=32GB

module load lammps/openmpi/intel/20210831

mpirun -np 32 --oversubscribe lmp -in in.aem
