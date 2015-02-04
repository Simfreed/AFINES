# README #

### What is this repository for? ###

* Coarse-grained model of actomyosin networks
<<<<<<< HEAD

### System Requirements ###
Minimally, this system requires gcc+11 and boost which you can load on midway via the commands:

    >>>module load gcc
    >>>module load boost

If you'd like to use my plotter, then you need Mathematica. 
Assuming you're accessing midway through ThinLinc, you can open mathematica through an interactive
computing node with the commands:

    >>>module load mathematica/9.0
    >>>mathematica

### QUICKSTART GUIDE ###

* If you don't already have a bin directory, create one with:
    
    >>> mkdir bin

* If you don't already have an executable, run the command: 
  
    >>> make network 

* you should now have an executable file called bin/nt. NOTE: you only need to recreate this file if you edit the source
  code.

* Create an output directory for your simulation (not necessarily named "output") as well as the "txt_stack" and "data"
  directories (necessarily named "txt_stack" and "data") e.g. with the commands:

    >>>mkdir output
    >>>mkdir output/txt_stack
    >>>mkdir output/data

* Run your simulation by executing your executable and specifying your output directory, e.g., 
    
    >>>./bin/nt --dir output

* See below for other simulation configuration options that you can set from the command line or from a configuration
  file

* Once your simulation has completed, the following files will have been generated:
    * output/txt_stack/rods.txt //the positions of every rod (actin) at every time step
    * output/txt_stack/links.txt //the positions of every link (these connect actin rods to make actin filaments) at
                                    every time step)
    * output/txt_stack/a_motors.txt //the positions of all active motors (e.g., myosin) at every time step
    * output/txt_stack/p_motors.txt //the positions of all passive motors (e.g., pmotors) at every time step
    * output/data/output.txt //some metadata about the simulation

    "position" in each case is a 4 tuple: (xcm, ycm, dx, dy)
    i.e., the center of mass position as well as the length in both dimensions of the particle. 
    Thus, the angle of the rod phi = ArcTan[dy/dx]
          it's length l = Sqrt[dy^2 + dx^2]
          the position of the pointed end = (xcm + l/2 Cos[phi], ycm + l/2 Sin[phi])
          the position fo the barbed  end = (xcm - l/2 Cos[phi], ycm - l/2 Sin[phi])
    
    The time associated with the positions is on it's own line before each list of positions within the file. 
    thus the structure of these files is:

    t = t1
    x1, y1, dx1, dy1
    .
    .
    .
    xn, yn, dxn, dyn
    t=t2
    x1, y1, dx1, dy1,
    .
    .
    .
    t=tn
    .
    .
    .
    xn, yn, dxn, dyn

* To plot the simulations, I wrote a mathematica script located at 

    analysis/simPlot.nb

 that can be opened within Mathematica. It was written using Mathematica 9.0, but I think will work with any version 8.0
 or higher. You're free to plot the simulations with your favorite plotter, in case you don't like Mathematica.  

### Configurable settings ###

Currently the following options for a simulation can be set upon execution, either from the command line, or within a
configuration file:

For example, to run a 500 second of simulation of 10 rigid actin filaments, an active motor density of 0.5 and a crosslinker density
of 0.05 you would enter the command:

    >>> ./bin/nt --tf 500 --npolymer 10 --a_motor_density 0.5 --p_motor_density 0.05

(this would write to the default output directory)

For an example configuration file, see below:

* xrange (double)  [default : 50] = size of cell in horizontal direction (um)
* yrange (double)  [default : 50] = size of cell in vertical direction (um)
         
* dt (double)  [default : 0.001] = length of individual timestep in seconds
* tfinal (double)  [default : 10] = length of simulation in seconds
* nframes (int)  [default : 1000] = number of frames of actin/link/motor positions printed to file (equally distant in
  time)
* nmsgs (int)  [default : 10000] = number of timesteps between printing simulation progress to stdout
        
* viscosity (double)  [default : 0.5] = Implicity viscosity to determine friction [um^2 / s]
* temperature,temp (double)  [default : 0.004] = Temp in kT [pN-um] that effects magnituded of Brownian component of simulation
* bnd_cnd,bc (string)  [default : "NONE"] = boundary conditions
         
* nmonomer (double)  [default : 1] = number of monomers per filament
* npolymer (double)  [default : 3] = number of polymers in the network
* actin_length (double)  [default : 10] = Length of a single actin monomer
* actin_pos_str (string)   [default : ""] = Starting positions of actin polymers, commas delimit coordinates; semicolons delimit positions
         
* a_motor_density (double)  [default : 0.001] = number of active motors / area
* p_motor_density (double)  [default : 0.001] = number of passive motors / area
* a_motor_pos_str (string)   [default : ""] = Starting positions of motors, commas delimit coordinates; semicolons delimit positions
         
* a_m_kon (double)  [default : 90.0] = active motor on rate
* a_m_koff (double)  [default : 1] = active motor off rate
* a_m_kend (double)  [default : 5] = active motor off rate at filament end
* a_motor_stiffness (double)  [default : 50] = active motor spring stiffness (pN/um)
         
* p_m_kon (double)  [default : 90] = passive motor on rate
* p_m_koff (double)  [default : 0] = passive motor off rate
* p_m_kend (double)  [default : 0] = passive motor off rate at filament end
* p_motor_stiffness (double)  [default : 50] = passive motor spring stiffness (pN/um)
         
* link_length (double)  [default : 0] = Length of links connecting monomers
* polymer_bending_modulus (double)  [default : 0.04] = Bending modulus of a filament
* fracture_force (double)  [default : 1000000] = pN-- filament breaking point
* bending_fracture_force (double)  [default : 1000000] = pN-- filament breaking point
* link_stretching_stiffness,ks (double)  [default : 10] = stiffness of link, pN/um
* use_linear_bending,linear (bool)  [default : true] =option to send spring type of bending springs
* shear_rate (double)  [default : 0] = shear rate in pN/(um-s)
         
* dir (string)  [default : "out/test"] = output directory
* seed (int)  [default : time(NULL] )= Random number generator seed

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
actin_length=10
a_motor_density=0.05
p_motor_density=0.5
actin_pos_str=0,0,0

/////Example ends at the line above////
=======
* 0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* run execute_simulation.py
* see execute_simulation.py
* g++, gnuplot, python
* How to run tests
* python execute_simulation.py
  pngs will be in appropriate amf_*_* folder
>>>>>>> c8a59e58ad872d99d9d09326e70f3522c5b6bd54

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Simon Freedman
<<<<<<< HEAD
* GCIS E126
=======
>>>>>>> c8a59e58ad872d99d9d09326e70f3522c5b6bd54
