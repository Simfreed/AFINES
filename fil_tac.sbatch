#!/bin/bash 
#SBATCH --job-name=fil_tac_sbatch 
#SBATCH --output=fil_tac_sbatch.out 
#SBATCH --error=fil_tac_sbatch.err
#SBATCH --time=36:00:00
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --cpus-per-task=1

./bin/afines --npolymer 1000 --p_motor_density 0.5 --xrange 5 --yrange 5 --nmonomer 5 --grid_factor 4 --tfinal 30 --dir ./out/exv_rc/ex19/

