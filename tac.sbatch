#!/bin/bash 
#SBATCH --job-name=fil_sbatch_2 
#SBATCH --output=fil_sbatch_2.out 
#SBATCH --error=fil_sbatch_2.err
#SBATCH --time=12:00:00
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --cpus-per-task=1

./bin/afines --npolymer 100 --p_motor_density 5 --xrange 5 --yrange 5 --tfinal 50 --link_length 0.5 --nmonomer 4 --grid_factor 4 --dir ./out/exv_rc/ex21/


