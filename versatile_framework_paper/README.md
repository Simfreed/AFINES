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
any simulation. I.e., to fullly run the simulation for figure 1D, 

 > mkdir slide
 > mkdir slide/txt_stack
 > mkdir slide/data
 > ./bin/afines -c versatile_framework_paper/config/slide.cfg --dir slide

To actually create any of the figure panels, see plots.nb



### Figure 1 ### 
* D

    > ./bin/afines -c versatile_framework_paper/config/slide.cfg

* E

     > ./bin/afines -c versatile_framework_paper/config/buckle.cfg


### Figure 2 ###
* A
 * Red
   
   > ./bin/afines -c versatile_framework_paper/config/contracting_short_lo_dens.cfg --npolymer [5*N]
   
   where N is 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000
 * Blue

    > ./bin/afines -c versatile_framework_paper/config/contracting_short.cfg --p_motor_density 0 --a_motor_density [0.001*N]

   where N is 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000
 * Black

    > ./bin/afines -c versatile_framework_paper/config/contracting_short.cfg --npolymer N --a_motor_density [0.1*N] --p_motor_density [0.1*N]

    where N is 1,2,5,10,20,50,100,200,500
* B 

    > ./bin/afines -c versatile_framework_paper/config/contracting_short.cfg --xrange [10*N] --yrange [10*N] --npolymer [5*N*N]

  where N is 1,2,3,4,5,6,7,8,9,10
* C

   > ./bin/afines -c versatile_framework_paper/config/contracting_short_hi_dens.cfg --grid_factor [0.01*N]

   where N is 9,10,20,30,40,50,60,70,80,90,100,200
### Figure 3 ###
* B
  
  > ./bin/afines -c versatile_framework_paper/config/L200_long.cfg --myseed [31415*N]
  
  where N is 1,2,...,100
* C
  
  > ./bin/afines -c versatile_framework_paper/config/one_hundred_filaments.cfg --myseed [31415*N] --actin_in versatile_framework_paper/config/one_hundred_filaments.txt
  
  where N is 1,2,...,100
### Figure 4 ###
* A
  
   ./bin/afines -c versatile_framework_paper/config/L200.cfg --dt 0.00001 --polymer_bending_modulus [0.01*N]
  
  where N is 1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100
* B
  
   ./bin/afines -c versatile_framework_paper/config/L200.cfg --dt 0.00001 --link_stretching_stiffness [0.01*N] 
  
  where N is 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000
* C
  
   ./bin/afines -c versatile_framework_paper/config/L200.cfg --dt 0.00001 --nmonomer N --link_length [200/N] --actin_length [100/N]
  
  where N is 10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000
### Figures 5 and S2 ###
 
  ./bin/afines -c versatile_framework_paper/config/shear.cfg --p_motor_stiffness [N*0.1]
 
 where N is 1,2,5,10,20,50,100,200,500,1000,2000,5000,10000
### Figures 6 and S3 ###
 * Green
  
   ./bin/afines -c versatile_framework_paper/config/motility.cfg --a_motor_density [N*0.1]
  
  where N is 0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80
 * Red
  
   ./bin/afines -c versatile_framework_paper/config/motility.cfg --nmonomer N
  
  where N is 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50
 * Blue
  
   ./bin/afines -c versatile_framework_paper/config/motility.cfg --a_m_kon [N*0.001]
  
  where N is 1,2,5,10,20,50,100,200,500,1000,2000
### Figures 7 and S4 ###
 
  ./bin/afines -c versatile_framework_paper/config/contracting.cfg
 
### Figure S1 ###
 
  ./bin/afines -c versatile_framework_paper/config/L200.cfg --nframes 20000 --tf 100
 
### Figure S2 ###
 
  ./bin/afines -c versatile_framework_paper/config/shear.cfg --n_bw_shear N --tf [0.00005*N]
 
 where N is 2,5,10,20,50,100,200,500,1000,2000,5000,10000
 NOTE: 0.00005 = g/dg * dt where g=total_strain=0.5, dg=0.001, dt=10^-7

