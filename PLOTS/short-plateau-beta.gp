#!/usr/bin/env gnuplot
#load "long-times.gp" #This calculates the right values of beta and plateau

reset


set term post enh c eps font "Times-Roman,40" size 6,6
set output "./FIGURES/short-plateau-beta.eps"

set multiplot layout 4,1
set lmargin 6
set rmargin 1
set tmargin 0
set bmargin 0
set xr[0.4:1.05]

#
# Short-time curvature
#
unset key
unset xlabel
set ylabel "{/Times-Italic a}_2" offset 1,0
set xtics auto format ""
set ytics 0,20 font ",28" offset 0.5,0
set label "(a)" at 0.9,10
plot [][0:50]"< awk '(NR%3==0)' ../OUTPUT/T*/N1080/shortNoise_NVT.txt"  using 1:5:6 w yerrorbars  ps 2 pt 7 lc rgb "dark-violet" title "{/Times-Italic K}({/Times-Italic t})",\
     "< awk '(NR%3==0)' ../OUTPUT/T*/N1080/shortDiag_NVT.txt"  using 1:5:6 w yerrorbars ps 2 pt 6 lc rgb "dark-green" t "{/Times-Italic C}_d({/Times-Italic t})"
#
# Plateau height
#
unset label
#set key center right font ",33" invert spacing 1.3 samplen .5 box
unset key
set ylabel "{/Times-Italic A}" offset 2,0
set ytics 0,0.1
unset xlabel
set label "(b)" at 0.9,.1
p [][0.04:0.35] "./SMALL-DATA/long-times_noise.txt" u 1:2:3 with yerrorbars ps 2 pt 7 lc rgb "dark-violet" title "{/Times-Italic K}({/Times-Italic t})",\
"./SMALL-DATA/long-times_diag.txt" u 1:2:3 with yerrorbars ps 2 pt 6 lc rgb "dark-green" title "{/Times-Italic C}_d({/Times-Italic t})"

# No errorbars:
# p [][0.04:0.35] "./SMALL-DATA/long-times_noise.txt" u 1:2 w p ps 2 pt 7 lc rgb "dark-violet" title "{/Times-Italic K}({/Times-Italic t})",\
#  "./SMALL-DATA/long-times_diag.txt" u 1:2 w p ps 2 pt 6 lc rgb "dark-green" title "{/Times-Italic C}_d({/Times-Italic t})"

#
# Stretching exponent
#
unset label
set ytics 0.5,0.1
set ylabel "{/Symbol-Oblique b}" offset 2,0
set ytics 0,0.2 font ",28" offset 0.5,0
set xtics font ",28" offset 0,0.5 format "%g"
#unset key
set key font ",30" spacing 1.2 samplen 1.5 box at 1.02,0.95
set label "(c)" at 0.9,.54
p [][.4:1.07] "./SMALL-DATA/long-times_noise.txt" u 1:6:7 w e ps 2 pt 7 lc rgb "dark-violet" title "{/Times-Italic K}({/Times-Italic t})",\
"./SMALL-DATA/long-times_diag.txt" u 1:6:7 w e ps 2 pt 6 lc rgb "dark-green" title "{/Times-Italic C}_d({/Times-Italic t})",\
"<awk '($1<0.6){print $0}' ../THERMALIZE/data/tauMCT_kmax28.txt" u 1:4:5 w e ps 2 pt 4 lc rgb "dark-orange" title "{/Times-Italic K}_{MCT}({/Times-Italic t})"

# No errorbars
# p [][.5:1.07] "./SMALL-DATA/long-times_noise.txt" u 1:6 w p ps 2 pt 7 lc rgb "dark-violet" title "{/Times-Italic K}({/Times-Italic t})",\
# "./SMALL-DATA/long-times_diag.txt" u 1:6 w p ps 2 pt 6 lc rgb "dark-green" title "{/Times-Italic C}_d({/Times-Italic t})"

set xlabel "{/Times-Italic T}" offset 0,5
set border 0
unset key
unset tics
unset ylabel
unset label
p [][-2:-1]100*x not

unset multiplot
set out
