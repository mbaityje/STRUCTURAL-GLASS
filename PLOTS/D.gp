#!/bin/env gnuplot
reset

#Tc from Kob-Andersen
Tc=0.435
#Tc from noise
A_noise=1
Tc_noise=0.43
expo_noise=1
fit A_noise/(x-Tc_noise)**expo_noise "../THERMALIZE/data/friction.txt" u ($1):($3) via A_noise,Tc_noise,expo_noise
#Tc from diagonal
A_diag=150000
Tc_diag=0.47
expo_diag=1.68
fit [:1] A_diag/(x-Tc_diag)**expo_diag "../THERMALIZE/data/friction.txt" u ($1):($4) via A_diag,Tc_diag,expo_diag
p "../THERMALIZE/data/friction.txt" u ($1):($4) w lp, A_diag/(x-Tc_diag)**expo_diag
p "../THERMALIZE/data/friction.txt" u ($1-Tc_diag):($4) w lp, A_diag/(x)**expo_diag




#For the moment I just use the average between the two
Tc=0.5*(Tc_noise+Tc_diag)


set term post enh c eps
set out "./FIGURES/friction.eps"
set title "Friction"
set logs
set ylabel "Friction"
set xlabel "Temperature"
p "../THERMALIZE/data/friction.txt" u ($1-Tc):($3) w lp title "Noise",\
"../THERMALIZE/data/friction.txt" u ($1-Tc):($4) w lp title "Diagonal",\
80/x**0.75 lw 2,\
2.0e5/x**1.65 lw 2
#comparare con Flenner,Szamel PRE 72 (2005)
#Sjogren(1980)


set out "./FIGURES/D.eps"
set title "Diffusion constant"
fmsd(x)=A_msd*(x-Tc_msd)**expo_msd
fdiag(x)=A_diag*(x-Tc_diag)**expo_diag
fnoise(x)=A_noise*(x-Tc_noise)**expo_noise
Tc_msd=0.435; expo_msd=1.75; A_msd=.05
Tc_diag=0.435; expo_diag=1.75; A_diag=.05
Tc_noise=0.435; expo_noise=1.75; A_noise=.05
fit [:.7] fmsd(x) "../THERMALIZE/data/D.txt" u 1:2:3 via A_msd,Tc_msd,expo_msd
fit [:.7] fdiag(x) "../THERMALIZE/data/friction.txt" u 1:(1./$4) via A_diag,Tc_diag,expo_diag
fit [:.7] fnoise(x) "../THERMALIZE/data/friction.txt" u 1:(1./$3) via A_noise,Tc_noise,expo_noise

Tc=(Tc_noise+Tc_diag+Tc_msd)/3
print "Tc_msd   =",Tc_msd,  "   expo_msd   =",expo_msd
print "Tc_diag  =",Tc_diag, "   expo_diag  =",expo_diag
print "Tc_noise =",Tc_noise,"   expo_noise =",expo_noise
print "<Tc>=",Tc


set logs
set ylabel "{/Times-Italic D}"
set xlabel "{/Times-Italic T-T}_c"
set ytics format "10^{%T}"
set key bottom right
p "../THERMALIZE/data/friction.txt" u ($1-Tc):(1./$3) w lp title "Noise",\
"../THERMALIZE/data/friction.txt" u ($1-Tc):(1./$4) w lp title "Diagonal",\
"../THERMALIZE/data/D.txt" u ($1-Tc):2:3 w errorl title "MSD"


set xlabel "{/Times-Italic T}"
unset logs
set logs y
p "../THERMALIZE/data/friction.txt" u ($1):(1./$3) w lp title "Noise",\
"../THERMALIZE/data/friction.txt" u ($1):(1./$4) w lp title "Diagonal",\
"../THERMALIZE/data/D.txt" u ($1):2:3 w errorl title "MSD"



