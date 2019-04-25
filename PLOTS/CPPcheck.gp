#!/usr/bin/env gnuplot
reset
set term post enh c eps font "Times-Roman,16" size 3,4

set logs x
set xlabel "{/Times-Italic t}"
set ylabel "{/Times-Italic C}^{P}(t)"
set tmargin 0
set rmargin 0.5
set lmargin 2
set bmargin 2

unset key
set xtics format "10^{%T}" nomirror font ",12" offset -.4,.4
set ytics nomirror font ",12" offset 0.7,0
unset xlabel; unset ylabel

set out "./FIGURES/CPPcheck.eps"
set multiplot layout 4,3
set label "{/Times-Italic T} = 5.0" at .05,4.5
p [:1][-1.5:5.5]"../OUTPUT/T5.0/N1080/CPPcheckT5.0_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label
set label "{/Times-Italic T} = 2.0" at .05,1.8
p [:1][-0.6:2.2]"../OUTPUT/T2.0/N1080/CPPcheckT2.0_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label
set label "{/Times-Italic T} = 1.0" at .05,0.9
p [:1][-0.3:1.1]"../OUTPUT/T1.0/N1080/CPPcheckT1.0_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label

set label "{/Times-Italic T} = 0.8" at .05,0.9*0.8
p [:1][-0.3*0.8:1.1*0.8]"../OUTPUT/T0.8/N1080/CPPcheckT0.8_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label
set label "{/Times-Italic T} = 0.7" at .05,0.9*0.7
p [:1][-0.3*0.7:1.1*0.7]"../OUTPUT/T0.7/N1080/CPPcheckT0.7_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label
set label "{/Times-Italic T} = 0.6" at .05,0.9*0.6
p [:1][-0.3*0.6:1.1*0.6]"../OUTPUT/T0.6/N1080/CPPcheckT0.6_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label

set label "{/Times-Italic T} = 0.55" at .05,0.9*0.55
p [:1][-0.3*0.55:1.1*0.55]"../OUTPUT/T0.55/N1080/CPPcheckT0.55_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label
set label "{/Times-Italic T} = 0.52" at .05,0.9*0.52
p [:1][-0.3*0.52:1.1*0.52]"../OUTPUT/T0.52/N1080/CPPcheckT0.52_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label
set label "{/Times-Italic T} = 0.49" at .05,0.9*0.49
p [:1][-0.3*0.49:1.1*0.49]"../OUTPUT/T0.49/N1080/CPPcheckT0.49_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label

set label "{/Times-Italic T} = 0.47" at .05,0.9*0.47
p [:1][-0.3*0.47:1.1*0.47]"../OUTPUT/T0.47/N1080/CPPcheckT0.47_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label
set label "{/Times-Italic T} = 0.46" at .05,0.9*0.46
p [:1][-0.3*0.46:1.1*0.46]"../OUTPUT/T0.46/N1080/CPPcheckT0.46_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label
set label "{/Times-Italic T} = 0.45" at .05,0.9*0.45
p [:1][-0.3*0.45:1.1*0.45]"../OUTPUT/T0.45/N1080/CPPcheckT0.45_NVT.txt" u 1:2 w lp lc rgb 'red' ps 0.2 t"Measured",\
"" u 1:3 w lp lc rgb 'blue' ps 0.2 t"From memory function", 0 lc 0 lt 0
unset label




unset multiplot
