### How to run simulations and analysis for ###
## Mechanical and Kinetics Factors Drive Sorting of Crosslinkers ##
### Simon Freedman, Cristian Suarez, Jonathan Winkelman, Gregory Voth, David Kovar, Aaron Dinner, Glen Hocky ###

This README file is a guide to reproducing the simulations and domain identification for the figures in the paper.  
Once these have been run, one can reproduce the specific plots, using the xlink_sorting_paper/plots.nb Mathematica Notebook. 

This version of AFINES used for this paper is algorithmically described in the Supplement of the paper. 
The code can be accessessed from the github AFINES repository via the command
```
git checkout growing_spacer_exv
```
Alternatively, from this xlink_sorting_paper directory
```
tar -xvf afines.tar
cd AFINES
```

Compile the program via the commands
```
make bundles
```
This will create the executable "AFINES/bin/bun"
Some possible compilation issues are addressed in the AFINES/README.md file.

Simulation trajectories are written to ${dir}/txt_stack with metadata in ${dir}/data where "dir" can be specified using
the --dir parameter. Further details for running simulations are available in the README.md file in the AFINES directory.

### Default parameters to be used below ###
```
cfg1 = config/no_grow.cfg
cfg2 = config/grow.cfg
fovx = 35
fovy = 20
nfixed = 9
nseg = 1
nrand = 166
l1 = 0.2
l2 = 0.3
kb0 = 0.04
kon = 2
koff = 0.05

x0      = -11.5
y0      = 0.375
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
s2dens  = 0.25
nrand   = ${s2dens}*${fovx}*${fovy} - ${nfixed}
for densrat in (0.25, 0.33, 0.5, 0.67, 0.8, 1, 1.2, 1.5, 2, 3, 4)
    for seed in (1,2,..., 20)*10^7
```
####  initialize system ####
```
s1dens  = ${densrat}*${s2dens}

dir  = out/f1/seed${seed}_s2dens${s2dens}
ddir = ${dir}/data
mkdir -p ${ddir}

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${l2} --nseg ${nseg}
python init/random_spacers.py ${ddir}/random.txt --fov ${fovx} ${fovy} --lm ${l2} --nm ${nrand} --seed ${seed}

cat ${ddir}/fixed.txt ${ddir}/random.txt > ${ddir}/s2in.txt
```
#### run simulation ####
```
bin/bun -c ${cfg1} --dir ${dir} --myseed ${seed} --spacer1_density ${s1dens} --spacer2_in ${ddir}/s2in.txt >> ${ddir}/out 2>> ${ddir}/err
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
for lendiff in (0,0.027,0.05,0.1,0.3,0.5,1)
    for seed in (1,2,..., 20)*10^7
```
####  initialize system with long domain ####
```
l2      = ${l1} + ${lendiff}
kb2     = ${kb0} * ${l2}
y0      = 0.25*(${l1} + ${l2}) 
apos    ="${x0},-${y0},${pi}:${x0},${y0},${pi}"

dir  = out/f2/afines/ldiff_s2init/seed${seed}/ldiff${lendiff}
ddir = ${dir}/data
mkdir -p ${ddir}

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${l2} --nseg ${nseg} 
python init/random_spacers.py ${ddir}/random.txt --fov ${fovx} ${fovy} --lm ${l2} --nm ${nrand} --seed ${seed}

cat ${ddir}/fixed.txt ${ddir}/random.txt > ${ddir}/s2in.txt
```
#### run simulation ####
```
bin/bun -c ${cfg1} --dir ${dir} --myseed ${seed} --spacer2_in ${ddir}/s2in.txt --actin_pos_str ${apos} \
        --spacer2_length ${l2} --s2bend ${kb2} >> ${ddir}/out 2>> ${ddir}/err
```
Convert to lattice and find domains as in Fig. 1

#### initialize system with short domain ####
same as ''initialize system with long domain'', except replace initialization with 

```
dir  = out/f2/afines/ldiff_s1init/seed${seed}/ldiff${lendiff}
ddir = ${dir}/data
mkdir -p ${ddir}

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${l1} --nseg ${nseg} 
python init/random_spacers.py ${ddir}/random.txt --fov ${fovx} ${fovy} --lm ${l1} --nm ${nrand} --seed ${seed}

cat ${ddir}/fixed.txt ${ddir}/random.txt > ${ddir}/s1in.txt
bin/bun -c ${cfg1} --dir ${dir} --myseed ${seed} --spacer1_in ${ddir}/s1in.txt --actin_pos_str ${apos} \
        --spacer2_length ${l2} --s2bend ${kb2} >> ${ddir}/out 2>> ${ddir}/err
```

Convert to lattice and find domains.

### Figure 2C (black) ###
```
mu = -0.008

for lendiff in (0,0.027,0.05,0.1,0.3,0.5,1)
    for seed in (1,2,..., 20)*10^7
```
####  initialize system ####
```
dir  = out/f2/lattice/ldiff/seed${seed}/ldiff${lendiff}
ddir = ${dir}/data
mkdir -p ${ddir}
```
#### run simulation ####
```
python lattice.py ${dir} --dL ${lendiff} --seed ${seed} --muA ${mu} --muB ${mu}
```
Find domains as in Fig. 1

### Figure 2D (red) ###
``` 
nrand = 0.25*${fovx}*${fovy} - ${nfixed}
for lkb in (0.01, 0.03, 0.068, 0.10, 0.3, 0.5) 
    for seed in (1,2,..., 20)*10^7
```
####  initialize system with long domain ####
```
dir  = out/f2/afines/lkb_s2init/seed${seed}/lkb${lkb}
ddir = ${dir}/data
mkdir -p ${ddir}

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${l2} 
python init/random_spacers.py ${ddir}/random.txt --fov ${fovx} ${fovy} --lm ${l2} --nm ${nrand} --seed ${seed}

cat ${ddir}/fixed.txt ${ddir}/random.txt > ${ddir}/s2in.txt
```
#### run simulation ####
```
bin/bun -c ${cfg1} --dir ${dir} --myseed ${seed} --spacer2_in ${ddir}/s2in.txt \
        --polymer_bending_modulus ${lkb} >> ${ddir}/out 2>> ${ddir}/err
```
Convert to lattice and find domains as in Fig. 1

####  initialize system with short domain ####
```
dir  = out/f2/afines/lkb_s1init/seed${seed}/lkb${lkb}
ddir = ${dir}/data
mkdir -p ${ddir}

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${l1} 
python init/random_spacers.py ${ddir}/random.txt --fov ${fovx} ${fovy} --lm ${l1} --nm ${nrand} --seed ${seed}

cat ${ddir}/fixed.txt ${ddir}/random.txt > ${ddir}/s1in.txt
bin/bun -c ${cfg1} --dir ${dir} --myseed ${seed} --spacer1_in ${ddir}/s1in.txt \
        --polymer_bending_modulus ${lkb} >> ${ddir}/out 2>> ${ddir}/err
```
Convert to lattice and find domains as in Fig. 1

### Figure 2D (black) ###
```
mu = -0.008
kT = 0.004
for lkb in (0.01, 0.03, 0.05, 0.068, 0.10, 0.15, 0.2, 0.25, 0.3, 0.5) 
    lp = ${lkb}/${kT}
    for seed in (1,2,..., 20)*10^7
```
####  initialize system ####
```
dir  = out/f2/lattice/lkb/seed${seed}/lkb${lkb}
ddir = ${dir}/data
mkdir -p ${ddir}
```
#### run simulation ####
```
python lattice.py ${dir} --Lp ${lp} --seed ${seed} --muA ${mu} --muB ${mu}
```
Find domains as in Fig. 1

### Figure 3 ###
Self contained in kmc.ipynb

### Figure 5A-B ###
```
s2len=0.227
for densrat in (0.1, 0.2, 0.5, 1, 1.5, 2, 3)
    for seed in (1,2,...,20)*10^7
        for initL (1, 2)
```
####  initialize system ####
```
s1dens = 0.25*${densrat}
dir  = out/f5/s2len${s2len}_initL${initL}/seed${seed}/s1dens${s1dens}
ddir = ${dir}/data
mkdir -p ${ddir}

lm = ${l1} if ${initL}==1 else ${l2}
inputParam = "spacer1_in" if ${initL}==1 else "spacer2_in"
nrandL = ${nrand} if ${initL}==2 else (${s1dens}*${fovx}*${fovy} - ${nfixed})

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${lm} 
python init/random_spacers.py ${ddir}/random.txt --fov ${fovx} ${fovy} --lm ${lm} --nm ${nrandL} --seed ${seed}

cat ${ddir}/fixed.txt ${ddir}/random.txt > ${ddir}/xl_in.txt
```
#### run simulation ####
```
bin/bun -c ${cfg2} --dir ${dir} --myseed ${seed} --${inputParam} ${ddir}/xl_in.txt \
        --spacer1_density ${s1dens} --lgrow ${lgrow} >> ${ddir}/out 2>> ${ddir}/err
```
Convert to lattice and find domains as in Fig. 1

Repeat for s2len=0.3

### Figure 5C-E ###
First two lgrows are only necessary for 5C
```
for lgrow in (0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1)
    for s1koff1 in (0.05, 0.1, 0.2)
        for seed in (1, 2,..., 20)*10^7
            for initL (1, 2)
```
####  initialize system ####
```
dir  = out/f5/s1koff1${s1koff1}_initL${initL}/seed${seed}/lgrow${lgrow}
ddir = ${dir}/data
mkdir -p ${ddir}

lm = ${l1} if ${initL}==1 else ${l2}
inputParam = "spacer1_in" if ${initL}==1 else "spacer2_in"

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${lm} 
python init/random_spacers.py ${ddir}/random.txt --fov ${fovx} ${fovy} --lm ${lm} --nm ${nrand} --seed ${seed}

cat ${ddir}/fixed.txt ${ddir}/random.txt > ${ddir}/xl_in.txt
```
#### run simulation ####
```
bin/bun -c ${cfg2} --dir ${dir} --myseed ${seed} --${inputParam} ${ddir}/xl_in.txt \
        --s1_koff1 ${s1koff1} --spacer1_density ${s1dens} --lgrow ${lgrow} >> ${ddir}/out 2>> ${ddir}/err
```
Convert to lattice and find domains as in Fig. 1
