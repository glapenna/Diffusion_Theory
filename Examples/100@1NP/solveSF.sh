#!/bin/bash
#SBATCH --account INF25_biophys_0
#SBATCH --partition=dcgp_usr_prod
##SBATCH --qos=dcgp_qos_dbg
#SBATCH --time=01:00:00
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-node=1
#SBATCH --mem=247000
##SBATCH --gres=gpu:1
#SBATCH --job-name viscosity

## si puo' arrivare fino a 2x56 task per nodo in dcgp

module purge
module load profile/base
module load openblas/0.3.24--gcc--12.2.0
module load openmpi/4.1.6--gcc--12.2.0                           
module load netlib-scalapack/2.2.0--openmpi--4.1.6--gcc--12.2.0  

export OMP_NUM_THREADS=1

time ../../bin/solveSF -n1 5900 -n2 0 -f fs_1.dump -o1 eval_1.out -o2 evec_m1.dump
