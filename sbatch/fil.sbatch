#!/bin/bash 
#SBATCH --job-name=fil_sbatch 
#SBATCH --output=fil_sbatch.out 
#SBATCH --error=fil_sbatch.err
#SBATCH --time=60:00:00
#SBATCH --qos=weare-dinner-long 
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --cpus-per-task=1

./bin/afines --npolymer 1000 --xrange 5 --yrange 5 --nmonomer 5 --grid_factor 4 --tfinal 25 --dir ./out/exv_rc/ex19/  
