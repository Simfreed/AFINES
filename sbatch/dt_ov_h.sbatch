#!/bin/bash
# Shell options


# set the number of nodes to 2
#SBATCH --job-name=ov_h_profile 
#SBATCH --output=/project/weare-dinner/ibunge/dt/ov_h/log/arrayJob_%A_%a.out
#SBATCH --error=/project/weare-dinner/ibunge/dt/ov_h/log/arrayJob_%A_%a.err
#SBATCH --array=0,1,2,3,4,5,6,7,8
#SBATCH --time=24:00:00
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --constraint=e5-2680v4
#SBATCH --cpus-per-task=1

###########################
# Begin work section      #
###########################

echo "Array Job; NOT REMAKING EXECUTABLE"
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

time=(0.001 0.00075 0.0005 0.00025 0.0001 0.000075 0.00005 0.000025 0.00001)

main_dr="/project/weare-dinner/ibunge/dt/ov_h/${time[$SLURM_ARRAY_TASK_ID]}"

data_dr="$main_dr/data"

cmd="--npolymer 100 --a_motor_density 0.0 --xrange 5 --yrange 5 --dt ${time[$SLURM_ARRAY_TASK_ID]} --grid_factor 2 --nmonomer 3 --link_length 0.5 --p_motor_density 1.0 --p_motor_length 0.10 --tfinal 50 --dir $main_dr"
echo $cmd
startt=`date +%s`
./bin/afines_ov $cmd
endt=`date +%s`
tott=$(echo $startt $endt | awk '{printf "%d", $2-$1}')
echo $tott > "$data_dr/time.txt"

mv "/project/weare-dinner/ibunge/dt/ov_h/log/arrayJob_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" $data_dr
mv "/project/weare-dinner/ibunge/dt/ov_h/log/arrayJob_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" $data_dr

