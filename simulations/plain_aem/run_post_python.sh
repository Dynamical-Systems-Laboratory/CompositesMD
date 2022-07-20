#!/bin/bash

#SBATCH --job-name=plainAEMPP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=48GB

module load python/intel/3.8.6

./run_python_part.sh
