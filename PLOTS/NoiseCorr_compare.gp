#!/usr/bin/env gnuplot
reset

set term post enh c eps font "Times-Roman,32"
########################################
# Show noise autocorrelation functions #
########################################

set output "./FIGURES/K_compare.eps"
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
plot[:][:]"../OUTPUT/T5.0/N1080/noisecorrJK_NVT_M3.txt"  using 1:2:3  lc rgb "red" w errorl title "Volterra+Laplace",\
"../OUTPUT/T5.0/N1080/noisecorrlinear_NVE.txt"  using 1:2 w lp title "Linear",\
"../OUTPUT/T5.0/N1080/noisecorrlinear-selfcon_NVE.txt"  using 1:2 w lp "blue" title "Selfcon"

set out
unset label


