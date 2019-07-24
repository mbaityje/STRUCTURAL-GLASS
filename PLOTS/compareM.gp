#!/usr/bin/env gnuplot

if (! defined(T)) T="5.0"

set title "T = ".T
set ylabel "{/Times-Italic t}"
set xlabel "{/Times-Italic C}^P({/Times-Italic t})"
set key bottom left
set logs x
set grid
p [][-.5:T]\
"../OUTPUT/T".T."/N1080/CPPcheck_NVT_M3.txt" u 1:2:3 w errorl title "{/Times-Italic M} = 3" lc rgb "#99775500",\
"../OUTPUT/T".T."//N1080/CPPcheck_NVT_M4.txt" u 1:2:3 w errorl title "{/Times-Italic M} = 4" lc rgb "#AA990055",\
"../OUTPUT/T".T."//N1080/CPPcheck_NVT_M5.txt" u 1:2:3 w errorl title "{/Times-Italic M} = 5" lc rgb "#CC009977",\
"../OUTPUT/T".T."//N1080/CPPcheck_NVT_M3.txt" u 1:4:5 w errorl notitle lc rgb "black" ps 0.1

#pause -1

reset
set title "T = ".T
set ylabel "{/Times-Italic t}"
set xlabel "{/Times-Italic K}({/Times-Italic t})"
set key bottom left
set logs x
set grid
p [][-.5:1]\
"../OUTPUT/T".T."/N1080/noisecorrJK_NVT_M3.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 3" lc rgb "blue" ps .5,\
"../OUTPUT/T".T."/N1080/noisecorrJK_NVT_M4.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 4" lc rgb "red" ps .5,\
"../OUTPUT/T".T."/N1080/noisecorrJK_NVT_M5.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 5" lc rgb "dark-green" ps .5

#pause -1

#
# TWO MULTIPLOTS
# 1) K(t) for all M at T=0.8
# 2) CPPcheck for all M at T=0.8
#



#
# K(t)
#

reset
set term post enh c eps font "Times-Roman,34"
set out "./FIGURES/compareM_T0.8.eps"
set multiplot
set tmargin at screen 0.985
set bmargin at screen 0.13
set rmargin at screen 0.98
set lmargin at screen 0.13
set ylabel "{/Times-Italic K}({/Times-Italic t})" offset 3,0
set xlabel "{/Times-Italic t}" offset 0,1.2
set key top right samplen 1 font ",28"
set tics font ",18"
set ytics offset 0.5,0
set xtics format "10^{%T}" offset 0,0.5
set logs x
left_tail=1.0
right_tail=1.4
top_tail=0.06
bottom_tail=0.04
set arrow nohead from left_tail,top_tail to left_tail,bottom_tail lw 0.4
set arrow nohead from right_tail,top_tail to right_tail,bottom_tail lw 0.4
set arrow nohead from left_tail,top_tail to right_tail,top_tail lw 0.4
set arrow nohead from left_tail,bottom_tail to right_tail,bottom_tail lw 0.4
#set arrow from 0.5*(left_tail+right_tail),top_tail to 2.2,0.28 nohead dt 2
set arrow from left_tail,top_tail to 0.65,0.28 nohead dt 2
set arrow from right_tail,top_tail to 7.45,0.28 nohead dt 2
left_bulk=0.075
right_bulk=0.4
top_bulk=0.23
bottom_bulk=0.085
set arrow nohead from left_bulk,top_bulk to left_bulk,bottom_bulk lw 0.4
set arrow nohead from right_bulk,top_bulk to right_bulk,bottom_bulk lw 0.4
set arrow nohead from left_bulk,top_bulk to right_bulk,top_bulk lw 0.4
set arrow nohead from left_bulk,bottom_bulk to right_bulk,bottom_bulk lw 0.4
set arrow nohead from left_bulk,bottom_bulk to 0.035,0.036 dt 2
set arrow nohead from left_bulk,top_bulk to 0.035,0.405 dt 2
p [:10][-.05:1]\
"../OUTPUT/T0.8/N1080/noisecorrJK_NVT_M3.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 3" lc rgb "blue" ps .5,\
"../OUTPUT/T0.8/N1080/noisecorrJK_NVT_M4.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 4" lc rgb "red" ps .5,\
"../OUTPUT/T0.8/N1080/noisecorrJK_NVT_M5.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 5" lc rgb "dark-green" ps .5,\
0 lc 0 dt 3 notitle
set tmargin at screen 0.6
set bmargin at screen 0.4
set rmargin at screen 0.95
set lmargin at screen 0.7
unset arrow
unset key
unset label; unset xlabel; unset ylabel
unset tics
p [left_tail:right_tail][bottom_tail:top_tail]\
"../OUTPUT/T0.8/N1080/noisecorrJK_NVT_M3.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 3" lc rgb "blue" ps .5,\
"../OUTPUT/T0.8/N1080/noisecorrJK_NVT_M4.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 4" lc rgb "red" ps .5,\
"../OUTPUT/T0.8/N1080/noisecorrJK_NVT_M5.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 5" lc rgb "dark-green" ps .5
set tmargin at screen 0.5
set bmargin at screen 0.2
set rmargin at screen 0.4
set lmargin at screen 0.2
p [left_bulk:right_bulk][bottom_bulk:top_bulk]\
"../OUTPUT/T0.8/N1080/noisecorrJK_NVT_M3.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 3" lc rgb "blue" ps .5,\
"../OUTPUT/T0.8/N1080/noisecorrJK_NVT_M4.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 4" lc rgb "red" ps .5,\
"../OUTPUT/T0.8/N1080/noisecorrJK_NVT_M5.txt" u 1:4:5 w errorl title "{/Times-Italic M} = 5" lc rgb "dark-green" ps .5

unset multiplot


#
# CPPcheck(t)
#

set term post enh c eps font "Times-Roman,34"
set out "./FIGURES/CPPcheckM_T0.8.eps"
set multiplot
set tmargin at screen 0.985
set bmargin at screen 0.13
set rmargin at screen 0.98
set lmargin at screen 0.13
set ylabel "{/Times-Italic C^P}({/Times-Italic t})" offset 3.2,0
set xlabel "{/Times-Italic t}" offset 0,1.3
set key top right samplen 1 font ",28"
set tics font ",18"
set ytics offset 0.5,0
set xtics format "10^{%T}" offset 0,0.5
set logs x
left_bulk=0.1
right_bulk=0.4
top_bulk=-0.01
bottom_bulk=-0.21
set arrow nohead from left_bulk,top_bulk to left_bulk,bottom_bulk lw 0.4
set arrow nohead from right_bulk,top_bulk to right_bulk,bottom_bulk lw 0.4
set arrow nohead from left_bulk,top_bulk to right_bulk,top_bulk lw 0.4
set arrow nohead from left_bulk,bottom_bulk to right_bulk,bottom_bulk lw 0.4
set arrow nohead from left_bulk,top_bulk to 0.0478,0.219 dt 2
set arrow nohead from left_bulk,bottom_bulk to 0.0478,-0.395 dt 2
p [:1][-.5:0.8]\
"../OUTPUT/T0.8/N1080/CPPcheck_NVT_M3.txt" u 1:2:3 w errorl title "{/Times-Italic M} = 3" lc rgb "#99775500",\
"../OUTPUT/T0.8/N1080/CPPcheck_NVT_M4.txt" u 1:2:3 w errorl title "{/Times-Italic M} = 4" lc rgb "#AA990055",\
"../OUTPUT/T0.8/N1080/CPPcheck_NVT_M5.txt" u 1:2:3 w errorl title "{/Times-Italic M} = 5" lc rgb "#CC009977",\
"../OUTPUT/T0.8/N1080/CPPcheck_NVT_M3.txt" u 1:4:5 w errorl title "Measured" lc rgb "black" ps 0.1,\
0 lc 0 dt 3 notitle
set tmargin at screen 0.6
set bmargin at screen 0.2
set rmargin at screen 0.55
set lmargin at screen 0.2
unset arrow
unset key
unset label; unset xlabel; unset ylabel
unset tics
set object 1 rectangle from screen 0.2,0.2 to screen 0.55,0.6
p [left_bulk:right_bulk][bottom_bulk:top_bulk]\
"../OUTPUT/T0.8/N1080/CPPcheck_NVT_M3.txt" u 1:2:3 w errorl title "{/Times-Italic M} = 3" lc rgb "#99775500",\
"../OUTPUT/T0.8//N1080/CPPcheck_NVT_M4.txt" u 1:2:3 w errorl title "{/Times-Italic M} = 4" lc rgb "#AA990055",\
"../OUTPUT/T0.8//N1080/CPPcheck_NVT_M5.txt" u 1:2:3 w errorl title "{/Times-Italic M} = 5" lc rgb "#CC009977",\
"../OUTPUT/T0.8//N1080/CPPcheck_NVT_M3.txt" u 1:4:5 w errorl title "Measured" lc rgb "black" ps 0.1,\
0 lc 0 dt 3 notitle
unset multiplot
