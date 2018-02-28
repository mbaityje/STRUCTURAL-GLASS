#!/usr/bin/gnuplot

set term png enh font "Times-Roman,24"

set tmargin 1
set bmargin 2
set rmargin 1
set lmargin 6

set ylabel "{/Times-Italic E}_{IS}" offset 3,0
set xlabel "Time step" offset 0,1
set tics font ",16"
set ytics offset .5,0
set xtics offset 0,.4
unset key

set out "./FIGURES/eIStraj_T049S0.png"
set ytics -400,4
p[:12599999][-400:-384] "../OUTPUT/T0.49/N65/S0/elistIS.txt" u 1:2 pt 8 ps .1 lc rgb '#8B0000' notitle

set out "./FIGURES/eIStraj_T049S1.png"
set ytics -402,4
p[:12799999][-402:-382] "../OUTPUT/T0.49/N65/S1/elistIS.txt" u 1:2 pt 8 ps .1 lc rgb '#8B0000' notitle

set out "./FIGURES/eIStraj_T049S1_zoom.png"
set rmargin 3
set ytics -400,1
set xtics 5255000,4000
p [5.255e6:5.263e6][-400:-397]"../OUTPUT/T0.49/N65/S1/elistIS.txt" u 1:2 w lp pt 8 ps .3 lc rgb '#8B0000' notitle
