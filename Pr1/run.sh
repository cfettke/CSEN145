#!/bin/bash
#
#SBATCH --job-name=pl-lab3
#SBATCH --partition=cpu
#CPU Cores = 2*omp_get_max_threads()
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#BATCH --output=pl-lab3-%j.out
#SBATCH --time=10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER@scu.edu

export OMP_NUM_THREADS=X
export OMP_PLACES=cores
export OMP_PROC_BIND=true

make
./run.out 4000 2000 4000 > out.txt
