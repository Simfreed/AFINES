#!/bin/bash 
#SBATCH --job-name=1000_ov_sbatch 
#SBATCH --output=1000_ov_sbatch.out 
#SBATCH --error=1000_ov_sbatch.err
#SBATCH --time=12:00:00 
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --cpus-per-task=1

./old_version --npolymer 1000 --xrange 5 --yrange 5 --tfinal 20 --grid_factor 4 --nmonomer 2 --link_length 0.25 --rmax 0.1 --dir ./out/exv_rc/ex30

./old_version --npolymer 1000 --xrange 5 --yrange 5 --tfinal 20 --grid_factor 4 --nmonomer 2 --link_length 0.25 --p_motor_density 0.5 --rmax 0.1 --dir ./out/exv_rc/ex31
