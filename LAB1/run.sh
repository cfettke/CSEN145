#!/bin/bash
#
#SBATCH --job-name=pl-lab1
#SBACTH --partition=cpu
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=1G
#SBATCH --nodes=1
#SBATCH --output=pl-lab1-%j.out
#SBATCH --time=10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfettkether@scu.edu

export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

rm mat_vec_mult
make mat_vec_mult
./mat_vec_mult
diff myresults.txt result.txt
