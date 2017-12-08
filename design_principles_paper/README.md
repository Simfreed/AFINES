# README for running the simulations in the paper                           #
# Design principles for selective self-assembly of active networks          #
# Simon L. Freedman, Shiladitya Banerjee, Glen M. Hocky, Aaron R. Dinner    #

This README file is a guide to reproducing the simulations that the figures in the paper analyze.  
IF you have not yet read the README file in the top directory, read that first                    

The code in design_principles_paper/afines.tar directory was used to produce the plots.         
This code is accessible from the github under the tag v1.0; i.e., if you downloaded AFINES from its github, you can access it from the top AFINES directory 
with the command
```
git checkout 1.0
```
or, to check it into a new branch
```
git checkout -b design_principles_version v1.0
```

Should things change and the current version of the code is incompatible, with these config files 
or gives weird results, try compiling and using that version of the code.                         

Simulation output will be written to the directories 
```
<dir>/txt_stack
<dir>/data
```
where <dir> can be specified using the "--dir" parameter and the default value is your working directory. 
<dir> must exist before the simulation runs. 


These instructions are for generating the raw trajectory data, and in some cases post-processing, but not to make the panels directly.
To actually create any of the figure panels, see design_principles_paper/mathematica/plots.nb.

### Executable Compilation ###
From inside the afines directory (i.e., if you enter "ls", you should minimally see, "src/", "prog/", "include/", "test/, and "makefile") compile the main dynamics package with 
```
make network
```
which should produce "bin/afines" (see README for issues) and the motor trajectory generator with
```
make motorwalk
```
which should produce "bin/motor_walk"

### Figures 1,2,3A-C,S1,S2 ### 
```
 bin/afines -c design_principles_paper/config/network.cfg --a_motor_density [M] --p_motor_densty [X] --myseed [S] --dir dens_dens/seed[S]_mdens[M]_xldens[X]
```
where: 
[S] is 1e7, 2e7, 3e7, 4e7, 5e7
[M] is 0.0, 0.01, 0.02 0.05, 0.1, 0.2, 0.3
[X] is 0.0, 0.1, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5
for a total of 5x7x8=280 simulations

### Figures 3D-F,S3,S5 ###
```
 bin/afines -c design_principles_paper/config/network.cfg --a_m_koff [M] --a_m_kend [M] --p_m_koff [X] --p_m_kend [X] --myseed [S] --dir koff_koff_len/seed[S]_len10_mkoff[M]_xlkoff[X]
```
where: 
[S] is 1e7, 2e7, 3e7, 4e7, 5e7
[M] is 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1
[X] is 0.01, 0.02, 0.05, 0.1, 0.2 ,0.5, 1
for a total of 5x7x7=245 simulations

### Figures 4A, left ###
```
 bin/afines -c design_principles_paper/config/network.cfg --a_m_koff [M] --a_m_kend [M] --p_m_koff 0.05 --p_m_kend 0.05 --myseed [S] 
 --nmonomer [NM] --npolymer [NP] --dir koff_koff_len/seed[S]_len[L]_mkoff[M]_xlkoff0.05
```
where: 
[S] is 1e7, 2e7, 3e7, 4e7, 5e7
[M] is 0.01, 0.02, 0.05, 0.1, 0.2 ,0.5, 1
[L] is 1, 2, 5, 8, 15
[NM] = [L] + 1
[NP] = 5000/[L]
for a total of 5x7x5=175 simulations

(Note: [L] = 10 was run in the 3D-F simulations.)
### Figures 4A, right ###
```
 bin/afines -c design_principles_paper/config/network.cfg --a_m_koff 0.2 --a_m_kend 0.2 --p_m_koff [X] --p_m_kend [X] --myseed [S] 
 --nmonomer [NM] --npolymer [NP] --dir koff_koff_len/seed[S]_len[L]_mkoff0.2_xlkoff[X]
```
where: 
[S] is 1e7, 2e7, 3e7, 4e7, 5e7
[X] is 0.01, 0.02, 0.1, 0.2 ,0.5, 1
[L] is 1, 2, 5, 8, 15
[NM] = [L] + 1
[NP] = 5000/[L]
for a total of 5x6x5=150 simulations

(Note: [L] = 10 should have been run in the 3D-F simulations and [X] = 0.05 should have been run in the 4A simulations.)
### Figure 4B, S4A, S6 ###
```
 bin/afines -c design_principles_paper/config/network.cfg --a_m_koff [M] --a_m_kend [M] --p_motor_density 0 --a_motor_density 0.3 --myseed [S] 
 --nmonomer [NM] --npolymer [NP] --dir mkoff_len/seed[S]_len[L]_mkoff[M]
```
where: 
[S] is 1e7, 2e7, 3e7, 4e7, 5e7
[M] is 0.01, 0.02, 0.05, 0.1, 0.2 ,0.5, 1
[L] is 1, 2, 5, 8, 10, 15
[NM] = [L] + 1
[NP] = 5000/[L]
for a total of 5x7x6=210 simulations

### Figures 4C, S4B, S7 ###
```
 bin/afines -c design_principles_paper/config/network.cfg --p_m_koff [X] --p_m_kend [X] --a_motor_density 0 --p_motor_density 1.5 --myseed [S] 
 --nmonomer [NM] --npolymer [NP] --dir xlkoff_len/seed[S]_len[L]_xlkoff[X]
```
where: 
[S] is 1e7, 2e7, 3e7, 4e7, 5e7
[X] is 0.01, 0.02, 0.05, 0.1, 0.2 ,0.5, 1
[L] is 1, 2, 5, 8, 10, 15
[NM] = [L] + 1
[NP] = 5000/[L]
for a total of 5x7x6=210 simulations

### Figures 5A, S8A,C, S9C ###
The following process should be done for each of the 20 directories [srcdir]:
```
dens_dens/seed[S]_mdens0.00_xldens0.0
dens_dens/seed[S]_mdens0.30_xldens0.0
dens_dens/seed[S]_mdens0.00_xldens1.5
dens_dens/seed[S]_mdens0.30_xldens1.5
```
where: 
[S] is 1e7, 2e7, 3e7, 4e7, 5e7

```
mkdir -p [srcdir]/txt_stack/restart_more_motors/in
tail -5000 [srcdir]/txt_stack/actins.txt > [srcdir]/restart_more_motors/actins.txt
mkdir -p [srcdir]/restart_more_motors/seed[Z]
bin/motor_walk --dir [srcdir]/restart_more_motors/seed[Z] --actin_in [srcdir]/restart_more_motors/in/actins.txt -c config/motors10.cfg --myseed [Z]
```
where:
[Z] is 1e6, 2e6, ..., 125e6

### Figures S8B, S9A,B ###
Assuming the first part is done, immediately above,
```
mkdir -p [srcdir]/restart_more_motors_kend10/seed[Z]
bin/motor_walk --dir [srcdir]/restart_more_motors_kend10/seed[Z] --actin_in [srcdir]/restart_more_motors/in/actins.txt -c config/motors10.cfg --myseed [Z] --a_m_kend 10
```
where:
[Z] is 1e6, 2e6, ..., 125e6


### Figure 5B ###
The following process should be done for each of the 20 directories [srcdir]:
```
dens_dens/seed[S]_mdens0.00_xldens0.0
dens_dens/seed[S]_mdens0.30_xldens0.0
dens_dens/seed[S]_mdens0.00_xldens1.5
dens_dens/seed[S]_mdens0.30_xldens1.5
```
where: 
[S] is 1e7, 2e7, 3e7, 4e7, 5e7

```
mkdir -p [srcdir]/restart_shear/in
tail -5000 [srcdir]/txt_stack/actins.txt > [srcdir]/restart_shear/in/actins.txt
```
if [srcdir] has mdens0.00
```
touch [srcdir]/restart_shear/amotors.txt
```
else
```
tail -750 [srcdir]/txt_stack/amotors.txt > [srcdir]/restart_shear/amotors.txt
```
if [srcdir] has xldens0.0
```
touch [srcdir]/restart_shear/pmotors.txt
```
else
```
tail -3750 [srcdir]/txt_stack/pmotors.txt > [srcdir]/restart_shear/pmotors.txt
```
Then,
```
bin/afines -c design_principles_paper/config/shear.cfg --dir [srcdir]/restart_shear \
--actin_in [srcdir]/restart_shear/in/actins.txt  --p_motor_in [srcdir]/restart_shear/in/pmotors.txt --a_motor_in [srcdir]/restart_shear/in/amotors.txt
 
```


### Order Parameter Calculations ###
If you have access to Mathematica and Python on command line, perform the following for all trajectories where ${dir} is
the trajectory output directory.
Set variables for the directory with the analysis folder scripts and prepare output directory:
```
afineshome="/home/simonfreedman/Code/cytomod"
math="/software/mathematica-10.2-x86_64/bin/MathematicaScript"
mathdir="${afineshome}/analysis/"
fovx="50"
fovy="50"
dr="0.05"
maxr="25"
ti=1
tf=400
nsamp=10000

mkdir ${dir}/analysis
```

#### g(r) ####
```
$math -script ${mathdir}/gr.m $dir $fovx $fovy $dr $maxr $ti $tf $nsamp
```
this should output the file "${dir}/analysis/gr_dr${dr}_maxr${maxr}_t${ti}-${tf}.dat"
To make most of the plots more efficiently, run
```
tail -50 ${dir}/analysis/gr_dr${dr}_maxr${maxr}_t1-400.dat > ${dir}/analysis/gr_dr${dr}_maxr${maxr}_t351-400.dat
```
#### g(r_barbed) ####
```
$math -script ${mathdir}/gr_barb.m $dir $fovx $fovy $dr $maxr $ti $tf $nsamp
```
this should output the file "${dir}/analysis/gr_barbed_dr${dr}_maxr${maxr}_t${ti}-${tf}.dat"
To make most of the plots more efficiently, run
```
tail -50 ${dir}/analysis/gr_barbed_dr${dr}_maxr${maxr}_t1-400.dat > ${dir}/analysis/gr_barbed_dr${dr}_maxr${maxr}_t351-400.dat
```
#### mesh size ####
This one's real slow.
```
voxel="0.25"
subdivs=10
emptythresh=1
dt=1
$math -script ${mathdir}/mesh_size.m $dir $ti $tf $voxel $subdivs $emptythresh $fovx $fovy $dt 
```
this should output the file "${dir}/analysis/meshsizes_1D_vx${voxel}_thresh${emptythresh}_t${ti}-${tf}.dat"
To make most of the plots more efficiently, run
```
tail -50 ${dir}/analysis/meshsizes_1D_vx${voxel}_thresh${emptythresh}_t1-400.dat > ${dir}/analysis/meshsizes_1D_vx${voxel}_thresh${emptythresh}_t351-400.dat
```
#### divergence ####
First calculate the velocity field with python
```
nb=10
dt=10
python ${mathdir}/rbf_vield.py $dir --nbins $nb --minpts 10 --dt $dt 
```
then run the mathematica script that processes it
```
thresh=10
dr=1
cmd="$math -script $mathdir/div_from_weights.m $main_dr $fovx $fovy $dr $dt $tf $nb $thresh"
```
#### mean squared displacement and angle distribution ####
For every s in 1..125, and every ${dir} used for figure 5:
```
${mdir}=${dir}/restart_more_motors_long/seed${s}e6
ti=5000
tf=30000
```
Unroll the trajectory
```
python ${mathdir}/periodic2definite.py ${mdir} amotors $fovx $fovy
```
Calculate angle distributions
```
$math -script ${mathdir}/vmot_ang_hists.m $mdir $ti $tf 0.001 Nothing 1 2 5 10 20 50 100 200 500 1000
```
Calculate msd
```
$math -script ${mathdir}/msd_tf_const.m $mdir $fovx $fovy $ti $tf
```
