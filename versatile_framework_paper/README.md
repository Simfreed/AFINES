# README for running the simulations in the paper                                                   #
# A versatile framework for simulating the dynamic mechanical structure of cytoskeletal networks    #
# Simon L. Freedman, Shiladitya Banerjee, Glen M. Hocky, Aaron R. Dinner                            #

This README file is a guide to reproducing the simulations that the figures in the paper analyze.  
It will not reproduce the analysis, only the trajectories that were analyzed.                     
To reproduce the analysis, see the versatile_framework_paper/plots.nb Mathematica Notebook       
IF you have not yet read the README file in the top directory, read that first                    

The code in versatile_framework_paper/afines.tar directory was used to produce the plots.         
Should things change and the current version of the code is incompatible, with these config files 
or gives weird results, try compiling and using that version of the code.                         

In all of these cases, you have to specify and create the output directory  just as you would for 
any simulation. I.e., to fully run the simulation for figure 1D, 

```
 mkdir slide
 mkdir slide/txt_stack
 mkdir slide/data
 ./bin/afines -c versatile_framework_paper/config/slide.cfg --dir slide
```

These instructions are for generating the raw data, and in some cases post-processing, but not make the panels directly.
To actually create any of the figure panels, see plots.nb.
For nearly every panel, there is a conversion from  "txt" to "dat" format required for making the plot in mathematica.
This can be done with an included python script using the following command:
```
python versatile_framework_paper/python/txt2dat.py [dir] [objects]
```
where [dir] is the top directory of the simulation output (e.g., "slide") in the above example, and [objects] is one of
"actins", "links", "amotors" or "pmotors".


### Figure 1 ### 
#### D ####
```
 ./bin/afines -c versatile_framework_paper/config/slide.cfg
```
#### E ####
```
 ./bin/afines -c versatile_framework_paper/config/buckle.cfg
```

### Figure 2 ###
#### A #### 
##### Red ##### 
```bash
./bin/afines -c versatile_framework_paper/config/contracting_short_r1.cfg --npolymer [N]
```   
where N is 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000
##### Blue ##### 
```bash
./bin/afines -c versatile_framework_paper/config/contracting_short_r1.cfg --a_motor_density [0.001*N]
```
where N is 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000
##### Black ##### 
```bash
./bin/afines -c versatile_framework_paper/config/contracting_short_r1.cfg --npolymer [N] --a_motor_density [0.01*N] --p_motor_density [0.01*N]
```
where N is 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000
#### B ####  
##### Blue ##### 
```
 ./bin/afines -c versatile_framework_paper/config/contracting_short.cfg --xrange [10*N] --yrange [10*N] --npolymer [5*N*N]
```
where N is 1,2,3,4,5,6,7,8,9,10
##### Red ##### 
```
 ./bin/afines -c versatile_framework_paper/config/contracting_short_hi_dens.cfg --grid_factor [0.01*N]
```
where N is 0,2,4,6,8,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450
### Figure 3 ###
#### B #### 
```
 ./bin/afines -c versatile_framework_paper/config/L20_x1.cfg
```
### Figure S2 ###
```
./bin/afines -c versatile_framework_paper/config/L20_x1.cfg --nframes 20000 --tf 100
```
#### C #### 
```
 ./bin/afines -c versatile_framework_paper/config/one_hundred_filaments.cfg --myseed [31415*N] --actin_in versatile_framework_paper/config/one_hundred_filaments_tf1000.txt
```
where N is 1,2,...,1000
Requires special post-processing, since this generates 250GB of data:
```
python python/periodic2definite.py [dir] actins 100 100
python python/process_scripts/txt2dat_ends.py [dir] actins_ext
```
where [dir] are the output directories from running the simulations.

### Figure 4 ###
#### A #### 
```
./bin/afines -c versatile_framework_paper/config/L20_x1.cfg --dt [D] --polymer_bending_modulus [0.001*N]
```
where N is 5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000
and 
D = 1e-4 if N <= 1000 ;
    1e-5 if 2000 <= N <= 50000 ; 
    1e-6 if N = 100000 ;

#### B #### 
```
./bin/afines -c versatile_framework_paper/config/L20_x1.cfg --dt [D] --link_stretching_stiffness [0.01*N] 
```
where N is 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000

and 
D = 1e-4 if N <= 1000 ;
    1e-5 if 2000 <= N <= 10000 ;
    1e-6 if N = 20000 ;

#### C #### 
```
./bin/afines -c versatile_framework_paper/config/L20_x1.cfg --dt [D] --link_length [0.01*N] --actin_length [0.005*N] --link_stretching_stiffness [100/N]
```
where N is 8,10,20,25,40,50,80,100,200,250,400,500
and 
D = 1e-4 if N >= 40 ; 
    1e-5 if 20 <= N <= 25 ; 
    1e-6 if N = 10 ; 

#### D #### 
```
./bin/afines -c versatile_framework_paper/config/L20_x1.cfg --dt [N*1e-6] 
```
where N is 1,2,5,10,20,50,100,200,500,1000,2000

### Figure 5 and S5 ###
#### Equilibration ####
```
./bin/afines -c versatile_framework_paper/config/shear_equil.cfg --p_motor_stiffness [N*0.1]
```
where N is 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000

#### Shear ####
```
tail -3168 {shear_equil_[N]}/txt_stack/pmotors.txt  in/pmotors.txt
tail -8000 {shear_equil_[N]}/txt_stack/actins.txt  in/actins.txt
./bin/afines -c versatile_framework_paper/config/shear.cfg --p_motor_stiffness [N*0.1] --actin_in in/actins.txt --p_motor_in in/pmotors/txt
```
where N is 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000
and {shear\_equil\_[N]} is the directory containing the corresponding equilibration

### Figure S3 ###
```
./bin/afines -c versatile_framework_paper/config/shear_intersected.cfg --n_bw_shear [N] --tf [N*0.00005] --actin_in in/actins.txt --p_motor_in in/pmotors.txt 
```
where N is 2,5,10,20,50,100,200,500,1000,2000,5000,10000
and in/actins.txt and in/pmotors.txt were generated from the N=10000 equilibration case above. 

#### NOTE: ####
* 0.00005 = g/dg * dt 
 * where g=total\_strain=0.5 
 * dg=0.001 
 * dt=10^-7

### Figure 6 ###
##### Green ##### 
```
./bin/afines -c versatile_framework_paper/config/motility.cfg --a_motor_density [N*0.01] --myseed [S] --dt [dt]
```
where N is 1,2,5,10,20,50,100,200,300,400,500,600,700,900,
S is 1e7, 2e7, 3e7, 4e7, 5e7,
and dt is 2e-5 if N>=300 and 5e-5 otherwise


##### Red ##### 
```
./bin/afines -c versatile_framework_paper/config/motility.cfg --nmonomer [N] --myseed [S]
```
where N is 2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,
S is 1e7, 2e7, 3e7, 4e7, 5e7


##### Blue ##### 
```
./bin/afines -c versatile_framework_paper/config/motility.cfg --a_m_kon [N*0.001] --dt [dt]
```
where N is 1,2,5,10,20,50,100,200,500,1000,2000,3000,4000,
S is 1e7, 2e7, 3e7, 4e7, 5e7,
and dt is 1e-5 if N is 3000, 2e-5 if N is 4000, and 5e-5 otherwise.

### Figures 7 ###
```
./bin/afines -c versatile_framework_paper/config/contracting.cfg --myseed [N]
```
where N is 1e7-20e7

### Figures S7 ###
```
./bin/afines -c versatile_framework_paper/config/contracting.cfg --myseed [N] --a_m_koff 1 --a_m_kend 10 --a_motor_density 1
```
where N is 1e7-20e7
#### 7D and S3D-G ####
Post processing to calculate the divergence
```
python versatile_framework_paper/python/div_by_vbin.py --nbins 10 --minpts 10 --padding 0.5 --dt 10
$MATHEMATICAHOME/bin/MathematicaScript -script versatile_framework_paper/div_from_weights.m [dir] 50 50 1 [dt] 399 10 10
```
where [dir] is an output directory of the simulation and [dt] is 1,2,5,10,20

### Figure S4 ###
#### A ####
```
/bin/sim versatile_framework_paper/config/actin_fiber20.cym
```

#### B ####
```
cp versatile_framework_paper/config/actin_fiber20.cym my_actin_fiber20.cym
sed -i "s/rigidity = 0.068/rigidity = [N*0.001]/g" my_actin_fiber20.cym
/bin/sim my_actin_fiber20.cym
```
where N is 5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000

#### C ####
```
cp versatile_framework_paper/config/actin_fiber20.cym my_actin_fiber20.cym
sed -i "s/segmentation = 1.00/segmentation = [N*0.01]/g" my_actin_fiber20.cym
/bin/sim my_actin_fiber20.cym
```
where N is 8,10,20,25,40,50,80,100,200,250,400,500

### Figure S6 ###
##### Green ##### 
```
cp config/motile.cym my_motile.cym
sed -i "s/10000 single particle/[N*25] single particle/g" my_motile.cym
sed -i "s/random_seed = 1000000/random_seed = [S]/g" my_motile.cym
/bin/sim my_motile.cym
```
where N is 0,1,2,5,10,20,50,100,200,300,400,500,600,700
and S is 1e6, 2e6, 3e6, 4e6, 5e6
##### Red ##### 
```
cp config/motile.cym my_motile.cym
sed -i "s/length = 15/length = N/g" my_motile.cym
sed -i "s/random_seed = 1000000/random_seed = [S]/g" my_motile.cym
/bin/sim my_motile.cym
```
where N is 1-9,11,13,15,17,19,21,23,25
and S is 1e6, 2e6, 3e6, 4e6, 5e6
##### Blue ##### 
```
cp config/motile.cym my_motile.cym
sed -i "s/ binding_rate = 1/ binding_rate = [N*0.001]/g" $cfg
sed -i "s/random_seed = 1000000/random_seed = [S]/g" $cfg
/bin/sim my_motile.cym
```
where N is 1,2,5,10,20,50,100,200,300,400,500,600,700,800,900,1000,2000,3000
and S is 8675309, 1e7, 2e7, 3e7, 4e7, 5e7
