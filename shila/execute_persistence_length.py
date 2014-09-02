import os
import subprocess
import sys

os.system('rm per')
use_fourier = raw_input("Use Fourier Modes to calculate persistence Length?\n")
if use_fourier.lower().startswith('y'):
    os.system('g++ persistence_length_fourier.cpp -o per')
    prefix = 'plf'
else:
    os.system('g++ persistence_length.cpp -o per')
    prefix = 'pl'

for i in range(1,20,2):
    for j in range(1,2,1):
        nmonomer = i
        link_stiffness = 10*j
        os.system('./per {0} {1}'.format(nmonomer, link_stiffness)) 

        main_dr = "persistence_length/{0}_kl{1}_nm{2}".format(prefix, link_stiffness, nmonomer) 
#        main_dr =raw_input('Enter a directory name for this simulation:\n')
#        main_dr = "amf_{0}_{1}".format(round(Actin_density,3), round(Motor_density,3))
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

