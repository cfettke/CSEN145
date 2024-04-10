#!/bin/bash
#
#SBATCH --job-name=LAB6
#SBATCH --partition=cpu
# CPU Cores = 2*omp_get_max_threads()
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=4
#SBATCH --output=pl-lab6-%j.out
#SBATCH --time=60:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfettkether@scu.edu

export OMP_NUM_THREADS=20
export OMP_PLACES=cores
export OMP_PROC_BIND=true

module load OpenMPI
mpirun ./lab6 ./example.clu > out.txt
