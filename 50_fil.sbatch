#!/bin/bash 
#SBATCH --job-name=50_fil_sbatch 
#SBATCH --output=50_fil_sbatch.out 
#SBATCH --error=50_fil_sbatch.err
#SBATCH --time=12:00:00
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --cpus-per-task=1

./bin/afines --npolymer 50 --p_motor_density 2 --xrange 5 --yrange 5 --tfinal 50 --link_length 0.5 --p_m_kon 10 --nmonomer 4 --grid_factor 4 --dir ./out/exv_rc/ex23/


