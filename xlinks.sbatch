#!/bin/bash 
#SBATCH --job-name=xlinks 
#SBATCH --output=xlinks.out 
#SBATCH --error=xlinks.err
#SBATCH --time=12:00:00
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --cpus-per-task=1

./bin/afines --npolymer 50 --p_motor_density 2.5 --p_m_kon 10 --xrange 5 --yrange 5 --tfinal 25 --link_length 0.5 --nmonomer 4 --grid_factor 4 --dir ./out/exv_rc/ex22/
