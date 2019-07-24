#!/usr/bin/env gnuplot
reset
set term post enh c eps font "Times-Roman,32"

set out "./FIGURES/friction.eps"
set logs
set ylabel "Friction"
set xlabel "Temperature"




#FITS FRICTION
A_noise=10; Tc_noise=0.435; expo_noise=2
fit [:.7] A_noise/(x-Tc_noise)**expo_noise "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" using ($1):($3):(0.1*$3) yerror via A_noise,Tc_noise,expo_noise
plot A_noise/(x-Tc_noise)**expo_noise, "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u ($1):($3)
replot A_noise/(x-Tc_noise)**expo_noise, "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoiseJK_NVT.txt" u ($1):($3):($4) w errorl

A_noise_long=10; Tc_noise_long=0.435; expo_noise_long=2
fit [:.7] A_noise_long/(x-Tc_noise_long)**expo_noise_long "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u ($1):($4):(0.1*$4) yerror via A_noise_long,Tc_noise_long,expo_noise_long
plot A_noise_long/(x-Tc_noise_long)**expo_noise_long, "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u ($1):($4)

A_diag=1500; Tc_diag=0.435; expo_diag=2
fit [:.7] A_diag/(x-Tc_diag)**expo_diag "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u ($1):($3):(0.1*$3) yerror via A_diag,Tc_diag,expo_diag
plot A_diag/(x-Tc_diag)**expo_diag, "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u ($1):($3)

A_diag_long=400; Tc_diag_long=0.435; expo_diag_long=2.25
fit [:.7] A_diag_long/(x-Tc_diag_long)**expo_diag_long "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u ($1):($4):(0.1*$4) yerror via A_diag_long,Tc_diag_long,expo_diag_long
plot A_diag_long/(x-Tc_diag_long)**expo_diag_long, "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u ($1):($4)


set output "./FIGURES/friction_all.eps"
load "~/.gnuplotting/gnuplot-palettes-master/paired.pal"
unset logs
set tmargin .5
set rmargin .5
set logs y
set xlabel "{/Times-Italic T}"
set ylabel "{/Symbol-Oblique z}"
set key font ",26" samplen 1
set ytics auto format "10^{%T}" font ",22" offset 0.5,0
set xtics auto font ",22" offset 0,0.4
p[:2.05][1:2e6] "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u 1:3 w lp lc rgb "dark-violet" ps 2 pt 7 t "{/Times-Italic K}({/Times-Italic t}) - Total",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u 1:4 w lp t "Long" lc rgb "dark-violet" ps 2 pt 5,\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u 1:5 w lp t "Short" lc rgb "dark-violet" ps 2 pt 9,\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u 1:3 w lp lc rgb "dark-green" ps 2 pt 6 t "{/Times-Italic C}_d({/Times-Italic t}) - Total",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u 1:4 w lp lc rgb "dark-green" ps 2 pt 4 t"Long",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u 1:5 w lp lc rgb "dark-green" ps 2 pt 8 t"Short"
#A_noise_long/(x-Tc_noise_long)**expo_noise_long not,\
#A_noise/(x-Tc_noise)**expo_noise ls 1 not,\
#A_diag_long/(x-Tc_diag_long)**expo_diag_long ls 4 not,\
#A_diag/(x-Tc_diag)**expo_diag ls 2 not,\
set output

set output "./FIGURES/friction_all_inset.eps"
set multiplot
load "~/.gnuplotting/gnuplot-palettes-master/paired.pal"
unset logs
set tmargin .2
set rmargin .2
set bmargin 2.3
set lmargin 5
set logs y
set xlabel "{/Times-Italic T}" offset 0,1
set ylabel "{/Symbol-Oblique z}" offset 1,0
set key font ",24" samplen 1 spacing 1 maxrows 3 at 1.45,16
set ytics auto format "10^{%T}" font ",22" offset 0.5,0
set xtics auto font ",22" offset 0,0.4 offset 0,0.5
p[:2.05][1:4e6] "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u 1:3 w lp lc rgb "dark-violet" ps 2 pt 7 t "{/Times-Italic K}({/Times-Italic t}) - Total",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u 1:4 w lp t "Long" lc rgb "dark-violet" ps 2 pt 5,\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u 1:5 w lp t "Short" lc rgb "dark-violet" ps 2 pt 9,\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u 1:3 w lp lc rgb "dark-green" ps 2 pt 6 t "{/Times-Italic C}_d({/Times-Italic t}) - Total",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u 1:4 w lp lc rgb "dark-green" ps 2 pt 4 t"Long",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u 1:5 w lp lc rgb "dark-green" ps 2 pt 8 t"Short"
set tmargin .5
set rmargin .7
set bmargin 10
set lmargin 19
set logs
set ylabel font ",22" offset 1.5,0
set xlabel "{/Times-Italic T-T}_{fit}" font ",22" offset 0,1.4
set tics font ",16"
unset key
p[:2.05][1:4e6] "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u ($1-Tc_noise):3 w lp lc rgb "dark-violet" ps 2 pt 7 t "{/Times-Italic K}({/Times-Italic t}) - Total",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u ($1-Tc_noise_long):4 w lp t "Long" lc rgb "dark-violet" ps 2 pt 5,\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u ($1-Tc_diag):3 w lp lc rgb "dark-green" ps 2 pt 6 t "{/Times-Italic C}_d({/Times-Italic t}) - Total",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u ($1-Tc_diag_long):4 w lp lc rgb "dark-green" ps 2 pt 4 t"Long",\
A_noise_long/(x)**expo_noise_long not dt 2 lc 0,\
A_noise/(x)**expo_noise not dt 1 lc 0,\
A_diag_long/(x)**expo_diag_long not dt 2 lc 0,\
A_diag/(x)**expo_diag not dt 1 lc 0

unset multiplot

set output

set print
pr "A_noise, Tc_noise, expo_noise"
pr A_noise, Tc_noise, expo_noise
pr "A_noise_long, Tc_noise_long, expo_noise_long"
pr A_noise_long, Tc_noise_long, expo_noise_long
pr "A_diag, Tc_diag, expo_diag"
pr A_diag, Tc_diag, expo_diag
pr "A_diag_long, Tc_diag_long, expo_diag_long"
pr A_diag_long, Tc_diag_long, expo_diag_long



set output
set title "Diffusion constant"
fmsd(x)=A_msd*(x-Tc_msd)**expo_msd
Tc_msd=0.435; expo_msd=1.75; A_msd=.05
fit [:.6] fmsd(x) "../THERMALIZE/data/D.txt" u 1:2:3 via A_msd,Tc_msd,expo_msd
plot fmsd(x), "../THERMALIZE/data/D.txt" u 1:2:3

fdiag(x)=A_diag*(x-Tc_diag)**expo_diag
Tc_diag=0.435; expo_diag=1.75; A_diag=.05
fit [:.6] fdiag(x) "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u 1:(1./$3):(0.1/$3) yerror via A_diag,Tc_diag,expo_diag
p fdiag(x), "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u 1:(1./$3)

fnoise(x)=A_noise*(x-Tc_noise)**expo_noise
Tc_noise=0.435; expo_noise=1.75; A_noise=.05
fit [:.6] fnoise(x) "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u 1:(1./$3):(0.1/$3) yerror via A_noise,Tc_noise,expo_noise
plot fnoise(x), "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u 1:(1./$3)

fdiag_long(x)=A_diag_long*(x-Tc_diag_long)**expo_diag_long
Tc_diag_long=0.435; expo_diag_long=1.75; A_diag_long=.05
fit [:.6] fdiag_long(x) "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u 1:(1./$4):(0.1*$4) yerror via A_diag_long,Tc_diag_long,expo_diag_long
plot fdiag_long(x), "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u 1:(1./$4)

fnoise_long(x)=A_noise_long*(x-Tc_noise_long)**expo_noise_long
Tc_noise_long=0.435; expo_noise_long=1.75; A_noise_long=.05
fit [:.6] fnoise_long(x) "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u 1:(1./$4):(0.1*$4) yerror via A_noise_long,Tc_noise_long,expo_noise_long
plot fnoise_long(x), "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u 1:(1./$4)

set print
Tc=(Tc_noise+Tc_diag+Tc_msd)/3
print "Tc_msd   =",Tc_msd,  "   expo_msd   =",expo_msd
print "Tc_diag  =",Tc_diag, "   expo_diag  =",expo_diag
print "Tc_diag_long  =",Tc_diag_long, "   expo_diag_long  =",expo_diag_long
print "Tc_noise =",Tc_noise,"   expo_noise =",expo_noise
print "Tc_noise_long =",Tc_noise_long,"   expo_noise_long =",expo_noise_long
print "<Tc>=",Tc



set out "./FIGURES/Dall_Tc.eps"
set logs
set ylabel "{/Times-Italic D}"
set xlabel "{/Times-Italic T-T}_c"
set ytics format "10^{%T}"
set key bottom right
p "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u ($1-Tc):(1./$3) w lp title "Noise",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u ($1-Tc):(1./$4) w lp title "Noise long",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u ($1-Tc):(1./$3) w lp title "Diagonal",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u ($1-Tc):(1./$4) w lp title "Diagonal long",\
"../THERMALIZE/data/D.txt" u ($1-Tc):2:3 w errorl title "MSD"


set out "./FIGURES/Dall.eps"
set xlabel "{/Times-Italic T}"
unset logs
set logs y
p "<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u ($1):(1./$3) w lp title "Noise",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionNoise_NVT.txt" u ($1):(1./$4) w lp title "Noise long",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u ($1):(1./$3) w lp title "Diagonal",\
"<awk '(NR%3==0)' ../OUTPUT/T*/N1080/frictionDiag_NVT.txt" u ($1):(1./$4) w lp title "Diagonal long",\
"../THERMALIZE/data/D.txt" u ($1):2:3 w errorl title "MSD"



