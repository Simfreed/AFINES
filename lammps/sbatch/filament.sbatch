#!/bin/sh

#SBATCH --job-name=lammps
#SBATCH --output=lammps-%j.out
#SBATCH --constraint=ib
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --partition=weare-dinner
#SBATCH --time=03:00:00

module load lammps-plumed
mpirun lmp_intelmpi<in.filament
