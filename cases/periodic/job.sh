#!/bin/bash
#SBATCH --job-name=periodic_u3r2
#SBATCH --output=periodic-u3r2.out
#SBATCH --error=periodic-u3r2.error
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --qos=np
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
module load prgenv/gnu
module load hpcx-openmpi
module load cmake/3.19.5
module load netcdf4/4.7.4
module load fftw/3.3.9
mpirun -np 16 ./dales4
