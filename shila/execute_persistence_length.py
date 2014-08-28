import os
import subprocess
import sys

for i in range(1,2,1):
    for j in range(0,3,1):
        nmonomer = 100 * i
        link_stiffness = pow(10,j)
        os.system('g++ persistence_length.cpp -o per')
        os.system('./per {0} {1}'.format(nmonomer, link_stiffness)) 
        main_dr = raw_input('Enter a directory name for this simulation:\n')
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

        
        #ffmpegstring=["ffmpeg","-r","1/0.0075","-i","time%01d.png","-c:v","libx264","-r","30","-pix_fmt","yuv420p",movie_file]
        #subprocess.call(ffmpegstring)
        #os.system('ffmpeg -r 1/0.01 -i time%01d.png -c:v libx264 -r 30 -pix_fmt yuv420p movie_file')                                                                                
