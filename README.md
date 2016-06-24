# README #

### What is this repository for? ###

* Coarse-grained model of actomyosin networks

### System Requirements ###
Minimally, this system requires gcc+11 and boost which you can load on midway via the commands:

    >module load gcc
    >module load boost

### QUICKSTART GUIDE ###

* If you don't already have a bin directory, create one with:
    
    >>> mkdir bin

* If you don't already have an executable, run the command: 
  
    >make [clean] [tar] network 
    
    [clean] will delete the old executable
    [tar] will generate the file tars/amxbd.tar.gz

* you should now have an executable file called bin/nt. NOTE: you only need to recreate this file if you edit the source
  code.

* Create an output directory for your simulation (not necessarily named "output") as well as the "txt_stack" and "data"
  directories (necessarily named "txt_stack" and "data") e.g. with the commands:

    >mkdir output
    >mkdir output/txt_stack
    >mkdir output/data

* Run your simulation in the specified output output directory, e.g., 
    
    >./bin/nt --dir output

* See below for other simulation configuration options that you can set from the command line or from a configuration
  file

* Once your simulation has completed, the following files will have been generated:
    * output/txt_stack/actins.txt //the trajectories of every actin bead
    * output/txt_stack/links.txt //the trajectories of every link 
    * output/txt_stack/amotors.txt //the trajectories of all active motors (e.g., myosin) at every time step
    * output/txt_stack/pmotors.txt //the trajectories of all passive motors (e.g., crosslinkers) at every time step
    * output/data/thermo.txt //the energies of actin filaments
    * output/data/output.txt //some metadata about the simulation

    All files are tab delimited
   txt_stack/actins.txt has the format

```
   x y r idx
        (x, y)  = bead position, 
        r  = bead radius 
        idx = index of filament that the bead is on
```

    txt_stack/links.txt has the format
    x y dx dy idx
        (x, y) = end position of link closer to the barbed end of filament 
        (x + dx, y + dy) = end position of link farther from barbed end 
        idx = index of filament the link is on

    txt_stack/amotors.txt and txt_stack/pmotors.txt have the format
    x y dx dy fidx0 fidx1 lidx0 lidx1
        (x, y) = position of head 0 of motor
        (x + dx, y + dy) = position of head 1 of motor
        fidx0 = index of filament that head 0 is attached to (-1 if not attached)
        fidx1 = index of filament that head 1 is attached to (-1 if not attached)
        lidx0 = index of link that head 0 is attached to (-1 if fidx0 = -1)
        lidx1 = index of link that head 1 is attached to (-1 if fidx1 = -1)

    data/thermo.txt has the format 
    KE PE TE idx
        KE = total v^2 of filament 
        PE = total potential energy of filament
        TE = KE + PE
        idx = filament index
    
    The time associated with the positions/energies is on it's own line before 
    each list of positions within the file. Thus the structure of actins.txt is:

    t = t1
    x1, y1, r1, idx1
    .
    .
    .
    xn, yn, rn, idxn
    t=t2
    x1, y1, r1, idx1,
    .
    .
    .
    t=tn
    .
    .
    .
    xn, yn, rn, idxn

    And similarly for other output files
### Configurable settings ###

Currently the following options for a simulation can be set upon execution, either from the command line, or within a
configuration file:

For example, to run a 500 second of simulation of 10 rigid actin filaments, an active motor density of 0.5 and a crosslinker density
of 0.05 you would enter the command:

    >>> ./bin/nt --tf 500 --npolymer 10 --a_motor_density 0.5 --p_motor_density 0.05

(this would write to the default output directory)

For an example configuration file, see below:

***********ENVIRONMENT**************

* xrange (double)  [default : 50um] = size of cell in horizontal direction
* yrange (double)  [50um] = size of cell in vertical direction 
* grid_factor (double) [1um^(-2)] = number of grid boxes 
* dt (double)  [0.0001s] = length of individual timestep
* tinit (double) [0s] = time that recording of simulation starts
* tfinal (double)  [10s] = length of simulation
* nframes (int)  [1000] = number of frames of actin/link/motor positions printed to file (equally distant in
  time)
* nmsgs (int)  [10000] = number of timesteps between printing simulation progress to stdout
* viscosity (double)  [0.001 mg/um*s]  = Dynamic viscosity
* temperature,temp (double)  [0.004 pN-um] = Temp in energy units 
* bnd_cnd,bc (string)  ["PERIODIC"] = boundary conditions
* dir (string) ["out/test"] = directory for output files
* myseed int (time(NULL)) = seed of random number generator
         
***********ACTIN PROPERTIES**********

* nmonomer (double)  [11] = number of beads per filament
* npolymer (double)  [3] = number of polymers in the network
* actin_length (double)  [0.5um] = Length of a single actin monomer
* actin_pos_str (string)   [""] = Starting positions of actin polymers, commas delimit coordinates; semicolons delimit positions
* link_length (double)  [0] = Length of links connecting monomers
* polymer_bending_modulus (double)  [0.04pn*um^2] = Bending modulus of a filament
* fracture_force (double)  [1000000pN] = filament breaking point
* bending_fracture_force (double)  [1000000pN] = filament breaking point
* link_stretching_stiffness,ks (double)  [1 pN/um] = stiffness of link
         

***********MOTOR PROPERTIES**********

* a_motor_density (double)  [0.05 um^(-2)] = number of active motors 
* a_motor_pos_str (string)   [""] = Starting positions of motors, commas delimit coordinates; semicolons delimit positions
* a_m_kon (double)  [100 s^(-1)] = active motor on rate
* a_m_koff (double)  [20 s^(-1)] = active motor off rate
* a_m_kend (double)  [20 s^(-1)] = active motor off rate at filament end
* a_motor_stiffness (double)  [10 (pN/um)] = active motor spring stiffness
* a_motor_length (double)  [0.4 um] = length of motor
* a_m_stall (double) [10pN] = stall force of motors
* a_m_break (double) [10pN] = rupture force of motors
* a_m_bind (double) [0.04pN*um] = binding energy
* a_motor_v (double) [1um/s] = velocity along filaments towards barbed end when attached
* motor_intersect_flag (boolean) [false] = if true, then motors are placed at filament intersections
* a_linkage_prob (double) [1] = probability that filaments are linked by a motor if motor_intersect_flag = true
* dead_head_flag (boolean) [false] = if true, then head [dead_head] of all motors remains stationary throughout sim
* dead_head (int) [0] = can be 0 or 1; head that remains stationary if dead_head_flag=true

***********XLINK PROPERTIES**********

* p_motor_density (double)  [0.05] = number of passive motors / um^2
* p_motor_pos_str (string)   [""] = Starting positions of xlinks, commas delimit coordinates; semicolons delimit positions
* p_m_kon (double)  [100 s^(-1)] = passive motor on rate
* p_m_koff (double)  [20 s^(-1)] = passive motor off rate
* p_m_kend (double)  [20 s^(-1)] = passive motor off rate at filament end
* p_motor_stiffness (double)  [50 s^(-1)] = xlink spring stiffness (pN/um)
* p_motor_length (double)  [0.4 s^(-1)] = length of xlink
* p_m_stall (double) [0pN] = stall force
* p_m_break (double) [10pN] = rupture force
* p_m_bind (double) [0.04pN*um] = binding energy
* link_intersect_flag (boolean) [false] = if true, then crosslinks are placed at filament intersections
* p_linkage_prob (double) [1] = probability that filaments are crosslinked if link_intersect_flag = true
* p_dead_head_flag (boolean) [false] = if true, then head [p_dead_head] of all xlinks remains stationary throughout sim
* p_dead_head (int) [0] = can be 0 or 1; head that remains stationary if p_dead_head_flag=true
* static_cl_flag (boolean) [false] = should be set to true if xlinks start off attached to filaments and never detach

*********SHEAR PROPERTIES**************

* strain_pct (double)  [0] = pre-strain (e.g., 0.5 means a strain of 0.5*xrange)
* time_of_strain (double) [0] = time of pre-strain

* d_strain_pct (double)  [0] = differential strain (e.g., 0.5 means a strain of 0.5*xrange)
* time_of_dstrain (double) [10000] = time when differential strain begins
* diff_strain_flag (boolean) [false] = flag to use if differential strain should be linear (in one direction)
* osc_strain_flag (boolean) [false] = flag to use if differential strain should be oscillatory (like Gardel, Science 2004)
* n_bw_shear (int) [10^8] = number of timesteps between subsequent differential strains 
* d_strain_freq (double) [1] = frequency of differential oscillatory strain

### Configuration file Example ###
The following is an example of a configuration file named example.cfg 
To use it together with a simulation, run the command
    
    >>>./bin/nt -c example.cfg

////Example begins below this line///

npolymer=500
nmonomer=1
dt=0.001
nframes=2000
tfinal=100
actin_length=0.5
a_motor_density=0.05
p_motor_density=0.5
actin_pos_str=0,0,0

/////Example ends at the line above////

### Contribution guidelines ###

* None yet, I should make some!
* Code review
* Other guidelines

### Who do I talk to? ###

* Simon Freedman
* GCIS E126