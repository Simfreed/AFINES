#!/bin/bash 
#SBATCH --job-name=tho 
#SBATCH --output=fil_sbatch_tho.out 
#SBATCH --error=fil_sbatch_tho.err
#SBATCH --time=12:00:00
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --cpus-per-task=1

./bin/afines --npolymer 200 --p_motor_density 5 --p_m_kon 10 --xrange 5 --yrange 5 --tfinal 50 --link_length 0.5 --nmonomer 4 --grid_factor 4 --dir ./out/exv_rc/ex19/


