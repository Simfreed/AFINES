#!/bin/bash
# Shell options


# set the number of nodes to 2
#SBATCH --job-name=nl_h_profile 
#SBATCH --output=/project/weare-dinner/ibunge/npol_1/nl_h/log/arrayJob_%A_%a.out
#SBATCH --error=/project/weare-dinner/ibunge/npol_1/nl_h/log/arrayJob_%A_%a.err
#SBATCH --array=1,2,5,10,20,50,100,200,500,1000
#SBATCH --time=70:00:00
#SBATCH --partition=weare-dinner
#SBATCH --account=weare-dinner
#SBATCH --constraint=e5-2680v4
#SBATCH --cpus-per-task=1

###########################
# Begin work section      #
###########################

echo "Array Job; NOT REMAKING EXECUTABLE"
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

np=$SLURM_ARRAY_TASK_ID

main_dr="/project/weare-dinner/ibunge/npol_1/nl_h/$np"

data_dr="$main_dr/data"

echo "deleting and recreating directory $main_dr"
rm -r $main_dr
mkdir -p $data_dr

cmd="--npolymer $np --a_motor_density 0.0 --xrange 5 --yrange 5 --dt 0.0001 --grid_factor 2 --nmonomer 2 --link_length 0.25 --p_motor_density 1.0 --p_motor_length 0.10 --tfinal 50 --dir $main_dr"
echo $cmd
startt=`date +%s`
./bin/afines $cmd
endt=`date +%s`
tott=$(echo $startt $endt | awk '{printf "%d", $2-$1}')
echo $tott > "$data_dr/time.txt"

mv "/project/weare-dinner/ibunge/npol_1/nl_h/log/arrayJob_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" $data_dr
mv "/project/weare-dinner/ibunge/npol_1/nl_h/log/arrayJob_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" $data_dr

