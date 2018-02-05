#!/bin/bash 
#SBATCH --job-name=db_sbatch 
#SBATCH --output=db_sbatch.out 
#SBATCH --error=db_sbatch.err
#SBATCH --time=24:00:00
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --constraint=e5-2680v4
#SBATCH --cpus-per-task=1

./bin/afines_debug --npolymer 1000 --xrange 5 --yrange 5 --nmonomer 3 --a_motor_density 0 --p_motor_density 0 --dt 0.00005 --grid_factor 4 --tfinal 25 --dir ./debug_test  
