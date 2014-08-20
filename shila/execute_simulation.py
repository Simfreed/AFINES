import os
import subprocess

for i in range(1,2,1):
    for j in range(1,2,1):
        npolymer = 1.
        nmonomer = 2.
        nmotor = 0 
        dx = 30.
        dy = 30.
        Mainpath = './'
        Output_Filename = "output"+str(i)+str(j)+".txt"
        Actin_length=3.
        Actin_density=npolymer * nmonomer/(dx*dy)
        Motor_length=0.5
        Motor_density=nmotor /(dx*dy)#i*0.1
        Motor_stiffness=50.
        Motor_speed=1.0
        Motor_kon=90.
        Motor_kend=5.
        Motor_koff=1.
        Actin_final="Actin_final_"+str(i)+str(j)+".txt"
        Myosin_final="Myosin_final_"+str(i)+str(j)+".txt"
        Viscosity=0.1+j*0.5
        movie_file = "movie"+str(i)+str(j)+".mp4"
        im_file = "im"+str(i)+str(j)+".png"
        os.system('g++ actin_myosin_flexible.cpp -o acto')
        os.system('./acto %s %s %f %f %f %f %f %f %f %f %f %s %s %f'%(Mainpath,\
                                                                   Output_Filename,\
                                                                   Actin_length,\
                                                                   Actin_density,\
                                                                   Motor_length,\
                                                                   Motor_density,\
                                                                   Motor_stiffness,\
                                                                   Motor_speed,\
                                                                   Motor_kon,\
                                                                   Motor_kend,\
                                                                   Motor_koff,\
                                                                   Actin_final,\
                                                                   Myosin_final,\
                                                                   Viscosity))
        main_dr = "amf_{0}_{1}".format(round(Actin_density,3), round(Motor_density,3))
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
