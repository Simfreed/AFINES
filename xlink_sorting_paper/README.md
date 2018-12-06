# How to reproduce plots for #
# Mechanical and Kinetics Factors Drive Sorting of Crosslinkers #
# Simon Freedman, Cristian Suarez, Jonathan Winkelman, Gregory Voth, David Kovar, Aaron Dinner, Glen Hocky #

This README file is a guide to reproducing the simulations that the figures in the paper analyze.  
It will not reproduce the analysis, only the trajectories that were analyzed.                     
To reproduce the analysis, see the xlink_sorting_paper/plots.nb Mathematica Notebook       
IF you have not yet read the README file in the top directory, read that first                    

The code in xlink_sorting_paper/afines.tar directory was used to produce the plots.         
Should things change and the current version of the code is incompatible, with these config files 
or gives weird results, try compiling and using that version of the code.                         

Simulation trajectories are written to ${dir}/txt_stack with metadata in ${dir}/data where "dir" can be specified using
the --dir parameter.

This version of AFINES used for this paper is algorithmically described in the Supplement of the paper. 
The code can be accessessed from the github via the command
```
git checkout growing_spacer_exv
```
Alternatively, from this xlink_sorting_paper directory
```
tar -xvf afines.tar
```
Compile via the commands
```
mkdir bin
make bundles
```
This will create the executable "bin/bun"

### Default parameters to be used below ###
```
cfg1 = config/no_grow.cfg
cfg2 = config/grow.cfg
fovx = 35
fovy = 20
nfixed = 9
l1 = 0.2
l2 = 0.3
kb0 = 0.04

x0      = -11.5
drx     = -0.1
dry     = 0

lmax = 15
lb = 0.037
ti = 1
tf = 2000
ncut = 4
rcut = 1

mmaex = /software/mathematica-10.2-x86_64/bin/MathematicaScript
```

### Figure 1 C-D ###
```
denstot = 0.5
for densrat in (0.25, 0.33, 0.5, 0.67, 0.8, 1, 1.2, 1.5, 2, 3, 4)
    for seed in (1,2,..., 20)*10^7
```
####  initialize system ####
```
s1dens  = ${densrat}*${denstot}/(1+${densrat})
s2dens  = ${denstot} - ${s1dens}
nrand   = ${s2dens}*${fovx}*${fovy} - ${nfixed}
x0      = -11.5
y0      = 0.25*(${l1} + ${l2}) 
apos    ="${x0},-${y0},${pi}:${x0},${y0},${pi}"

dir  = out/f1/seed${seed}_s2dens${s2dens}
ddir = ${dir}/data
mkdir -p ${ddir}

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${l2} 
python init/random_spacers.py ${ddir}/random.txt --fov ${fovx} ${fovy} --lm ${l2} --nm ${nrand} --seed ${seed}

cat ${ddir}/fixed.txt ${ddir}/random.txt > ${ddir}/s2in.txt
```
#### run simulation ####
```
bin/bun -c ${cfg1} --dir ${dir} --myseed ${seed} --spacer1_density ${s1dens} --spacer2_in ${ddir}/s2in.txt --actin_pos_str ${apos} >> ${ddir}/out 2>> ${ddir}/err
```
#### convert to lattice ####
```
mkdir -p ${dir}/analysis
${mmaex} -script process/afines2lattice.m ${dir} ${ti} ${tf} ${lb} ${lmax} 1 0 ${fovx} ${fovy}
```
#### find domains ####
```
${mmaex} -script process/domain_length.m ${dir} ${ti} ${tf} ${ncut} ${rcut} ${lb} ${lmax}
```

### Figure 2C (red) ###
``` 
nrand = 0.25*${fovx}*${fovy} - ${nfixed}
for lendiff in (0.000,0.010,0.027,0.05,0.1,0.2,0.5,1)
    for seed in (1,2,..., 20)*10^7
```
####  initialize system ####
```
l2      = ${l1} + ${lendiff}
kb2     = ${kb0} * ${l2}
y0      = 0.25*(${l1} + ${l2}) 
apos    ="${x0},-${y0},${pi}:${x0},${y0},${pi}"

dir  = out/f2/seed${seed}_ldiff${lendiff}
ddir = ${dir}/data
mkdir -p ${ddir}

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${l2} 
python init/random_spacers.py ${ddir}/random.txt --fov ${fovx} ${fovy} --lm ${l2} --nm ${nrand} --seed ${seed}

cat ${ddir}/fixed.txt ${ddir}/random.txt > ${ddir}/s2in.txt
```
#### run simulation ####
```
bin/bun -c ${cfg1} --dir ${dir} --myseed ${seed} --spacer2_in ${ddir}/s2in.txt --actin_pos_str ${apos} \
        --spacer2_length ${l2} --s2bend ${kb2} >> ${ddir}/out 2>> ${ddir}/err
```
Convert to lattice and find domains as in Fig. 1
### Figure 2D (red) ###
``` 
nrand = 0.25*${fovx}*${fovy} - ${nfixed}
for lkb in (0.01, 0.03, 0.05, 0.068, 0.10, 0.15, 0.2, 0.25, 0.3, 0.5) 
    for seed in (1,2,..., 20)*10^7
```
####  initialize system ####
```
apos    ="${x0},-${y0},${pi}:${x0},${y0},${pi}"

dir  = out/f2/seed${seed}_lkb${lkb}
ddir = ${dir}/data
mkdir -p ${ddir}

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${l2} 
python init/random_spacers.py ${ddir}/random.txt --fov ${fovx} ${fovy} --lm ${l2} --nm ${nrand} --seed ${seed}

cat ${ddir}/fixed.txt ${ddir}/random.txt > ${ddir}/s2in.txt
```
#### run simulation ####
```
bin/bun -c ${cfg1} --dir ${dir} --myseed ${seed} --spacer2_in ${ddir}/s2in.txt --actin_pos_str ${apos} \
        --polymer_bending_modulus ${lkb} >> ${ddir}/out 2>> ${ddir}/err
```
Convert to lattice and find domains as in Fig. 1

