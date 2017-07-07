#!/bin/bash 
#SBATCH --job-name=100_np_sbatch 
#SBATCH --output=100_np_sbatch.out 
#SBATCH --error=100_np_sbatch.err
#SBATCH --time=12:00:00 
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --cpus-per-task=1

./bin/afines --npolymer 100 --xrange 5 --yrange 5 --tfinal 50 -- nmonomer 5 --grid_factor 4 --dir ./out/exv_rc/ex26 

