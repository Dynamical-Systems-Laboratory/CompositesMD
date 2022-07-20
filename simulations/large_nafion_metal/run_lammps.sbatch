#!/bin/bash

#SBATCH --job-name=nf2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --time=5-00:00:00
#SBATCH --mem=8GB

module load lammps/openmpi/intel/20210831

mpirun -np 48 --oversubscribe lmp -in in.nafion_interface
