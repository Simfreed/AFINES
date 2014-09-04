import os
import subprocess
import sys

os.system('rm acto')
os.system('g++ simulation.cpp -o acto')

for i in range(1,2,1):
    for j in range(1,2,1):
        npolymer = 500 
        nmonomer = 10.
        nmotor = 1 
        actin_density=0.65 #npolymer * nmonomer/(dx*dy)
        motor_density=0.25 #nmotor /(dx*dy)#i*0.1
#        viscosity=0.1+j*0.5
        
        os.system('./acto {0} {1}'.format(actin_density, motor_density))
        
#        main_dr = raw_input('Enter a directory name for this simulation:\n')
        main_dr = "full_simulation/rhoA{0}_rhoM{1}".format(round(actin_density,2), round(motor_density,2))
        png_dr = main_dr + os.sep + "png_stack"
        txt_dr = main_dr + os.sep + "txt_files"

        print "gnuplotting"
        os.system('gnuplot "plotfiles.gp"')
        print "organizing output"
        os.system('rm -r ' + main_dr)
        os.system("mkdir {0}".format(main_dr))
        os.system('mkdir {0}'.format(png_dr))
        os.system('mkdir {0}'.format(txt_dr))
        os.system('mv *.png {0}'.format(png_dr))
        os.system('mv *.txt {0}'.format(txt_dr))

