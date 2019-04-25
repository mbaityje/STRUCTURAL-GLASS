#!/usr/bin/env gnuplot
reset

set term post enh c eps font "Times-Roman,32"
set out "./FIGURES/Dcheck.eps"
set multiplot
set lmargin 6
set tmargin .7
set rmargin 1
set bmargin 3
set xlabel "{/Times-Italic T}"
set ylabel "{/Times-Italic D}" offset 3,0
set key font ",28" at 5,0.05
set ytics 0,0.05 font ",24" nomirror offset .5,0
set xtics font ",24" nomirror offset 0,.5
set xr[0:5]
set yr[-0.005:]
p "../THERMALIZE/data/D.txt" u 1:2:3 w errorl t "from MSD",\
"<awk '$1==0' ../OUTPUT/T*/N1080/CPPcheckT*_NVT.txt| sort -nk7" u 7:4:5 w errorl t "from {/Times-Italic C}^{ {/Times-Italic P}}({/Times-Italic t})",\
"<awk '$1==0' ../OUTPUT/T*/N1080/CPPcheckT*_NVT.txt| sort -nk7" u 7:6 w lp t "from {/Times-Italic K}({/Times-Italic t})",\
0 lc 0 lt 0 not

set tmargin .7
set lmargin 10
set rmargin 15
set bmargin 10
unset key
unset ylabel
set ytics 0,0.01
set xtics 0.4,0.2
set xr[.4:1]
set yr[:0.018]
rep

unset multiplot
set out
reset
