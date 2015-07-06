#!/bin/sh

#SBATCH --job-name=lammps_filament
#SBATCH --output=/home/simonfreedman/scratch-midway/cytomod-lammps/log/langevin_%A_fric_%a.out
#SBATCH --error=/home/simonfreedman/scratch-midway/cytomod-lammps/log/langevin_%A_fric_%a.err
#SBATCH --constraint=ib
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# SBATCH --partition=weare-dinner
#SBATCH --account=pi-dinner
#SBATCH --time=10:00:00
#SBATCH --array=10,20,50,100,200,500,1000,2000,5000

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

bond_length=1 #um
natom=$(echo "$ideal_filament_length $bond_length" | awk '{printf "%d", $1/$2}') # = floor(i_f_l / b_l)
filament_length=$(echo "($natom - 1)*$bond_length" | bc)
mass=$(echo "scale=7; $mass_density * $filament_length / $natom" | bc)

timestep="0.0001"
filename="init_cfgs/filament_bl${bond_length}.txt"
angle_stiffness=$(echo "scale=7; $kT*$Lp/( 2*$bond_length )" | bc)  #"2e-5"
bond_stiffness="1e-2"
thermo_freq=10000
dump_freq=10000

vseed=8
lseed=8
pseed=8

friction=$(echo "$SLURM_ARRAY_TASK_ID" | awk '{printf "%.2f",$1/1000}')
dump_file="/home/simonfreedman/scratch-midway/cytomod-lammps/out/langevin_fric${friction}.atom"
iters=100000000
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

