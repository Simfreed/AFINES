do for [i=0:9999] {
   str_num=sprintf('%d',i)
   set palette maxcolors 3
   set palette defined ( 0.25 "yellow", 0.5 "green", 1 "blue")
   set term pngcairo enhanced 
   set style arrow 10 heads back nofilled linetype 3 linecolor palette linewidth 2.000 size screen 0.008,90.000,90.000 
   
   #set title "L_a=3.5, d_a=0.6, d_m=0.1, k_end=5.0"
   
   set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb"black" behind
   unset colorbox
   set xrange[-25:25]
   set yrange[-25:25]
   set output 'time'.str_num.'.png'
#   stats 'mfile'.str_num.'.txt' using 5 nooutput
   set cbrange[0.99:1.01]
   plot 'afile'.str_num.'.txt' u 1:2:3:4 with vectors filled head lw 2 linecolor rgb "red" notitle,\
   'lfile'.str_num.'.txt' u 1:2:3:4:5 notitle with vectors arrowstyle 10 ,\
   'amfile'.str_num.'.txt' u 1:2:3:4:5 notitle with vectors arrowstyle 10 ,\
   'pmfile'.str_num.'.txt' u 1:2:3:4:5 notitle with vectors arrowstyle 10
}
