#!/bin/sh

#SBATCH --job-name=lammps_filament
#SBATCH --output=/home/simonfreedman/scratch-midway/cytomod-lammps/log/filament_%A_T_%a.out
#SBATCH --error=/home/simonfreedman/scratch-midway/cytomod-lammps/log/filament_%A_T_%a.err
#SBATCH --constraint=ib
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# SBATCH --partition=weare-dinner
#SBATCH --account=pi-dinner
#SBATCH --time=10:00:00
#SBATCH --array=0-10,10-100:10,100-1000:100,1000-10000:1000

###########################
# Begin work section      #
###########################
echo "Array Job; NOT REMAKING EXECUTABLE"
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

temp=$SLURM_ARRAY_TASK_ID
kT300="0.000004" #pg-um^2/us^2
kT=$(echo "$temp $kT300" | awk '{printf "%.9f",$1*$2/300}') #pg-um^2/us^2

filament_length=200 #um
mass_density="0.000026" #pg/um
Lp="10" #um

bond_length=1
natom=200
mass=$(echo "$filament_length $mass_density $natom" | awk '{printf "%.6f",$1*$2/$3}')

timestep="0.0001"
filename="init_cfgs/filament_bl${bond_length}.txt"
angle_stiffness=$(echo "scale=7; $kT*$Lp/( 2*$bond_length )" | bc)  #"2e-5"
bond_stiffness="1e-2"
thermo_freq=10000
dump_freq=10000

dump_file="/home/simonfreedman/scratch-midway/cytomod-lammps/out/dump_filament_T${temp}.atom"
iters=100000000
exc="lmp_serial"
script="in/in.filament"

# Make the input file:
pythoncmd="python 2dfilament.py $bond_length"
echo $pythoncmd
$pythoncmd

echo $exc \
    -v timestep $timestep \
    -v filename $filename \
    -v temp $temp \
    -v my_mass $mass \
    -v angle_stiffness $angle_stiffness \
    -v bond_stiffness $bond_stiffness \
    -v bond_length $bond_length \
    -v thermo_freq $thermo_freq \
    -v dump_freq $dump_freq \
    -v dump_file $dump_file \
    -v iters $iters < $script

./$exc \
    -v timestep $timestep \
    -v filename $filename \
    -v temp $temp \
    -v my_mass $mass \
    -v angle_stiffness $angle_stiffness \
    -v bond_stiffness $bond_stiffness \
    -v bond_length $bond_length \
    -v thermo_freq $thermo_freq \
    -v dump_freq $dump_freq \
    -v dump_file $dump_file \
    -v iters $iters < $script

