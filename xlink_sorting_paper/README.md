### How to run simulations and analysis for ###
## Mechanical and Kinetics Factors Drive Sorting of Crosslinkers ##
### Simon Freedman, Cristian Suarez, Jonathan Winkelman, Gregory Voth, David Kovar, Aaron Dinner, Glen Hocky ###

This README file is a guide to reproducing the simulations and domain identification for the figures in the paper.  
Once these have been run, one can reproduce the specific plots, using the xlink_sorting_paper/plots.nb Mathematica Notebook. 

This version of AFINES used for this paper is algorithmically described in the Supplement of the paper. 
The code can be accessessed from the github via the command
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

dir  = out/f1/seed${seed}_s2dens${s2dens}
ddir = ${dir}/data
mkdir -p ${ddir}

python init/row_of_spacers.py ${ddir}/fixed.txt --r0 ${x0} ${y0} --dr ${drx} ${dry} --nm ${nfixed} --lm ${l2} 
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
for lendiff in (0.000,0.010,0.027,0.05,0.1,0.2,0.5,1)
    for seed in (1,2,..., 20)*10^7
```
####  initialize system ####
```
l2      = ${l1} + ${lendiff}
kb2     = ${kb0} * ${l2}
y0      = 0.25*(${l1} + ${l2}) 
apos    ="${x0},-${y0},${pi}:${x0},${y0},${pi}"

dir  = out/f2/afines/seed${seed}_ldiff${lendiff}
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

### Figure 2C (black) ###
```
mu = -0.008

for lendiff in (0.000,0.010,0.027,0.05,0.1,0.2,0.5,1)
    for seed in (1,2,..., 20)*10^7
```
####  initialize system ####
```
dir  = out/f2/lattice/seed${seed}_ldiff${lendiff}
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
for lkb in (0.01, 0.03, 0.05, 0.068, 0.10, 0.15, 0.2, 0.25, 0.3, 0.5) 
    for seed in (1,2,..., 20)*10^7
```
####  initialize system ####
```
dir  = out/f2/afines/seed${seed}_lkb${lkb}
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
dir  = out/f2/lattice/seed${seed}_lkb${lkb}
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
for s1dens in (0.25, 0.375, 0.5)
    for seed in (1,2,...,5)*10^7
        for initL (1, 2)
```
####  initialize system ####
```

dir  = out/f5/seed${seed}_initL${initL}_s1dens${s1dens}
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

### Figure 5D ###
```
for lgrow in (0.01, 0.1)
    for s1dens in (0.25, 0.375, 0.5)
        for s1koff1 in (0.05, 0.1, 0.2, 0.5)
            for seed in (1, 2, 3, 4, 5)*10^7
                for initL (1, 2)
```
####  initialize system ####
```
dir  = out/f5/seed${seed}_initL${initL}_lgrow${lgrow}_s1dens${s1dens}_s1koff1${s1koff1}
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
        --s1_koff1 ${s1koff1} --spacer1_density ${s1dens} --lgrow ${lgrow} >> ${ddir}/out 2>> ${ddir}/err
```
Convert to lattice and find domains as in Fig. 1

### Figure 5E ###
```
for lgrow in (0.01, 0.1)
    for s1dens in (0.25, 0.375, 0.5)
        for koffrat in (2, 4)
            for seed in (1, 2, 3, 4, 5)*10^7
                for initL (1, 2)
```
####  initialize system ####
```
s2koff = ${koff}/${koffrat}
s2kon  = ${kon}/${koffrat}
dir  = out/f5/seed${seed}_initL${initL}_lgrow${lgrow}_s1dens${s1dens}_koffrat${koffrat}
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
        --s2_kon ${s2kon} --s2_koff1 ${s2koff} --s2_koff2 ${s2koff} --s2_kend ${s2koff} \ 
        --spacer1_density ${s1dens} --lgrow ${lgrow} >> ${ddir}/out 2>> ${ddir}/err
```
Convert to lattice and find domains as in Fig. 1





