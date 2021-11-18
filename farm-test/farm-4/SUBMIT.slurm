#!/bin/bash --login

#SBATCH --job-name=farm-4
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

# Replace [budget code] below with your project code (e.g. t01)
#SBATCH --account=ecseab04
#SBATCH --partition=standard
#SBATCH --qos=standard

# Setup the batch environment
module load epcc-job-env

CASTEP=~/castep/obj/linux_x86_64_cray-XT--mpi/castep.mpi

# srun launches the parallel program based on the SBATCH options
srun --distribution=block:block --hint=nomultithread $CASTEP initial