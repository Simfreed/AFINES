#!/bin/bash 
#SBATCH --job-name=1000_np_sbatch 
#SBATCH --output=1000_np_sbatch.out 
#SBATCH --error=1000_np_sbatch.err
#SBATCH --time=12:00:00 
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --cpus-per-task=1

./bin/afines --npolymer 1000 --xrange 5 --yrange 5 --tfinal 50 --grid_factor 4 --nmonomer 5 --dt 0.001 --dir ./out/exv_rc/ex27
