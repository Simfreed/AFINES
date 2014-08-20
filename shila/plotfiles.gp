do for [i=0:7999] {
   str_num=sprintf('%d',i)  #Write integer to string
   str_num2=sprintf('%d',i) 
   set term png enhanced 
   set arrow 1 from -500, -100, 0 to 500, -100, 0 head back filled linetype 1 linecolor rgb "dark-violet"  linewidth 2.000 size screen 0.025,30.000,45.000
   set arrow 2 from -500, -110, 0 to 500, -110, 0 head back nofilled linetype 3 linecolor rgb "#56b4e9"  linewidth 2.000 size screen 0.030,15.000,90.000
   set arrow 3 from -500, -120, 0 to 500, -120, 0 head back filled linetype 1 linecolor rgb "dark-violet"  linewidth 2.000 size screen 0.030,15.000,45.000
   set arrow 4 from -500, -130, 0 to 500, -130, 0 head back filled linetype 3 linecolor rgb "#56b4e9"  linewidth 2.000 size screen 0.030,15.000,90.000
   set arrow 5 from -500, -140, 0 to 500, -140, 0 heads back filled linetype 1 linecolor rgb "dark-violet"  linewidth 2.000 size screen 0.030,15.000,135.000
   set arrow 6 from -500, -150, 0 to 500, -150, 0 head back empty linetype 3 linecolor rgb "#56b4e9"  linewidth 2.000 size screen 0.030,15.000,135.000
   set arrow 7 from -500, -160, 0 to 500, -160, 0 nohead back nofilled linetype 1 linecolor rgb "dark-violet"  linewidth 2.000
   set arrow 8 from -500, -170, 0 to 500, -170, 0 heads back nofilled linetype 3 linecolor rgb "#56b4e9"  linewidth 2.000 size screen 0.008,90.000,90.000
   set style line 1  linetype 1 linecolor rgb "dark-violet"  linewidth 2.000 pointtype 1 pointsize default pointinterval 0
   set style line 2  linetype 3 linecolor rgb "#56b4e9"  linewidth 2.000 pointtype 3 pointsize default pointinterval 0
   set style arrow 1 head back filled linetype 1 linecolor rgb "dark-violet"  linewidth 2.000 size screen 0.025,30.000,45.000
   set style arrow 2 head back nofilled linetype 3 linecolor rgb "#56b4e9"  linewidth 2.000 size screen 0.030,15.000,90.000
   set style arrow 3 head back filled linetype 1 linecolor rgb "dark-violet"  linewidth 2.000 size screen 0.030,15.000,45.000
   set style arrow 4 head back filled linetype 3 linecolor rgb "#56b4e9"  linewidth 2.000 size screen 0.030,15.000,90.000
   set style arrow 5 heads back filled linetype 1 linecolor rgb "dark-violet"  linewidth 2.000 size screen 0.030,15.000,135.000
   set style arrow 6 head back empty linetype 3 linecolor rgb "#56b4e9"  linewidth 2.000 size screen 0.030,15.000,135.000
   set style arrow 7 nohead back nofilled linetype 1 linecolor rgb "dark-violet"  linewidth 2.000
   set style arrow 8 heads back nofilled linetype 3 linecolor rgb "green"  linewidth 2.000 size screen 0.008,90.000,90.000
   #set title "L_a=3.5, d_a=0.6, d_m=0.1, k_end=5.0"
   set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb"black" behind
   set xrange[-25:25]
   set yrange[-25:25]
   set output 'time'.str_num2.'.png'
   plot 'afile'.str_num.'.txt' u 1:2:3:4 with vectors filled head lw 2 linecolor rgb "red" notitle,\
   'mfile'.str_num.'.txt' u 1:2:3:4 notitle with vectors arrowstyle 8 
}