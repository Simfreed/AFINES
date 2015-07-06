#!/bin/sh

#SBATCH --job-name=lammps_filament
#SBATCH --output=/home/simonfreedman/scratch-midway/cytomod-lammps/log/langevin_%A_bl_%a.out
#SBATCH --error=/home/simonfreedman/scratch-midway/cytomod-lammps/log/langevin_%A_bl_%a.err
#SBATCH --constraint=ib
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# SBATCH --partition=weare-dinner
#SBATCH --account=pi-dinner
#SBATCH --time=10:00:00
#SBATCH --array=1-40

###########################
# Begin work section      #
###########################
echo "Array Job; NOT REMAKING EXECUTABLE"
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

temp=300
kT300="0.000004" #pg-um^2/us^3
kT=$(echo "$temp $kT300" | awk '{printf "%.9f",$1*$2/300}') #pg-um^2/us^2

ideal_filament_length=200 #um
mass_density="0.000026" #pg/um
Lp="10" #um

bond_length=$(echo "$SLURM_ARRAY_TASK_ID" | awk '{printf "%.2f",$1/4}') #um
natom=$(echo "$ideal_filament_length $bond_length" | awk '{printf "%d", $1/$2}') # = floor(i_f_l / b_l)
filament_length=$(echo "($natom - 1)*$bond_length" | bc)
mass=$(echo "scale=7; $mass_density * $filament_length / $natom" | bc)

timestep="0.0001"
filename="init_cfgs/filament_bl${bond_length}.txt"
angle_stiffness=$(echo "scale=7; $kT*$Lp/( 2*$bond_length )" | bc)  #"2e-5"
bond_stiffness="1e-2"
thermo_freq=10000
dump_freq=1000

vseed=81
lseed=81
pseed=81

mu="0.0000005" #(pg/(um-us) Dynamic Viscosity 
fric_factor="10" #amount to multiply friction by
pi="3.1415926535897"
# Damping parameter in lammps:
# damp = mass/(4 pi mu R)
# R = bond_length/2 ==> damp = mass / (2 pi mu bond_length)

friction=$(echo "$fric_factor $mass $pi $mu $bond_length" | awk '{printf "%.2f", $1 * $2/(2 * $3 * $4 * $5)}')
dump_file="/home/simonfreedman/scratch-midway/cytomod-lammps/out/langevin_bl${bond_length}_fric${friction}.atom"
iters=10000000
exc="lmp_serial"
script="in/in.filament"

# Make the input file:
pythoncmd="python 2dfilament.py $bond_length $pseed"
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
    -v fric $friction \
    -v iters $iters \
    -v vseed $vseed \
    -v lseed $lseed \
    < $script

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
    -v fric $friction \
    -v vseed $vseed \
    -v lseed $lseed \
    -v iters $iters < $script

