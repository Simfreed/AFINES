#!/bin/bash 
#SBATCH --job-name=1000_fil_sbatch 
#SBATCH --output=1000_fil_sbatch.out 
#SBATCH --error=1000_fil_sbatch.err
#SBATCH --time=12:00:00 
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --cpus-per-task=1

./bin/afines --npolymer 1000 --p_motor_density 2.5 --xrange 5 --yrange 5 --tfinal 20 --p_m_kon 10 --grid_factor 4 --nmonomer 5 --dt 0.001 --dir ./out/exv_rc/ex25
