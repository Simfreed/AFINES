#!/bin/bash
# Shell options


# set the number of nodes to 2
#SBATCH --job-name=ov_profile 
#SBATCH --output=/project/weare-dinner/ibunge/ov/log/arrayJob_%A_%a.out
#SBATCH --error=/project/weare-dinner/ibunge/ov/log/arrayJob_%A_%a.err
#SBATCH --array=1,2,5,10,20,50,100,200,500,1000
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --account=pi-dinner
# SBATCH --qos=weare-dinner-long
# SBATCH --time=80:00:00
#SBATCH --cpus-per-task=1

###########################
# Begin work section      #
###########################

echo "Array Job; NOT REMAKING EXECUTABLE"
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

np=$SLURM_ARRAY_TASK_ID

main_dr="/project/weare-dinner/ibunge/ov/$np"

data_dr="$main_dr/data"

cmd="--npolymer $np --xrange 5 --yrange 5 --grid_factor 4 --nmonomer 3 --link_length 0.5 --p_motor_density 0.5 --p_m_kon 10 --tfinal 50 --dir $main_dr"
echo $cmd
startt=`date +%s`
./old_version $cmd
endt=`date +%s`
tott=$(echo $startt $endt | awk '{printf "%d", $2-$1}')
echo $tott > "$data_dr/time.txt"

mv "/project/weare-dinner/ibunge/ov/log/arrayJob_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" $data_dr
mv "/project/weare-dinner/ibunge/ov/log/arrayJob_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" $data_dr

