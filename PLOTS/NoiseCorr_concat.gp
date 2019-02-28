#!/usr/bin/env gnuplot
reset

# I plot all the temperatures, in case at some point I want to see them altogether
set logs x
plot[][-.1:1.1]"../OUTPUT/T5.0/N1080/noisecorrJK_NVT_M3.txt"  using 1:4:5 w errorl title "{/Times-Italic T} = 5.0  " ls 1 pt 1,\
"../OUTPUT/T2.0/N1080/noisecorrJK_NVT_M3.txt"  using 1:4:5 w errorl title "{/Times-Italic T} = 2.0  " ls 2 pt 1,\
"../OUTPUT/T1.0/N1080/noisecorrJK_NVT_M3.txt"  using 1:4:5 w errorl title "{/Times-Italic T} = 1.0  " ls 3 pt 1,\
"../OUTPUT/T0.8/N1080/noisecorrJK_NVT_M3.txt"  using 1:4:5 w errorl title "{/Times-Italic T} = 0.8  " ls 4 pt 1,\
"../OUTPUT/T0.7/N1080/noisecorrJK_NVT_M3.txt"  using 1:4:5 w errorl title "{/Times-Italic T} = 0.7  " ls 5 pt 1,\
"../OUTPUT/T0.6/N1080/noisecorrJK_NVT_M3.txt"  using 1:4:5 w errorl title "{/Times-Italic T} = 0.6  " ls 6 pt 1,\
"../OUTPUT/T0.55/N1080/noisecorrJK_NVT_M3.txt" using 1:4:5 w errorl title "{/Times-Italic T} = 0.55" ls 7 pt 1,\
"../OUTPUT/T0.52/N1080/noisecorrJK_NVT_M3.txt" using 1:4:5 w errorl title "{/Times-Italic T} = 0.52" ls 8 pt 1,\
"../OUTPUT/T0.49/N1080/noisecorrJK_NVT_M3.txt" using 1:4:5 w errorl title "{/Times-Italic T} = 0.49" ls 9 pt 1,\
"../OUTPUT/T0.47/N1080/noisecorrJK_NVT_M3.txt" using 1:4:5 w errorl title "{/Times-Italic T} = 0.47" ls 10 pt 1,\
"../OUTPUT/T0.46/N1080/noisecorrJK_NVT_M3.txt" using 1:4:5 w errorl title "{/Times-Italic T} = 0.46" ls 11 pt 1,\
"../OUTPUT/T0.45/N1080/noisecorrJK_NVT_M3.txt" using 1:4:5 w errorl title "{/Times-Italic T} = 0.45" ls 11 pt 1,\
0 ls 0 notitle



set term post enh c eps font "Times-Roman,32"

########################################
# Show noise autocorrelation functions #
########################################

set output "./FIGURES/K_concat.eps"
set lmargin 6
set tmargin .5
set rmargin .5
set logs x
set xlabel "time [LJ units]"
set ylabel "Memory function" offset 2.5,0
set key invert font ",27" top center
set tics font ",24"
set xtics format "10^{%T}"
load '~/.gnuplotting/gnuplot-palettes-master/rdbu.pal'
plot[:100][-50:500]"../OUTPUT/T0.52/N1080/noisecorrJK_NVT_M3.txt"  using 1:6:7 lc rgb "red" w errorl title "Volterra",\
"../OUTPUT/T0.52/N1080/noisecorrJK_NVT_M3.txt"  using 1:10:11 lc rgb "blue" w errorl title "Laplace"

set out
unset label


