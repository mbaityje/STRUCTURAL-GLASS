#!/usr/bin/env gnuplot

########################################
# Show noise autocorrelation functions #
########################################
reset

set title "Normalized Noise Autocorrelations for various T"
set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation"
set key invert
load '~/.gnuplotting/gnuplot-palettes-master/moreland.pal'
load '~/.gnuplotting/gnuplot-palettes-master/rdbu.pal'
plot   "../OUTPUT/T5.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "T = 5.0" ls 1
replot "../OUTPUT/T2.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "T = 2.0" ls 2
replot "../OUTPUT/T1.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "T = 1.0" ls 3
replot "../OUTPUT/T0.6/N1080/noisecorr_NVE.txt" using 1:3 w lp title "T = 0.6" ls 4
replot 0 ls 0 notitle

set title "Unnormalized Noise Autocorrelations for various T (with factor 1/T)"
plot   "../OUTPUT/T5.0/N1080/noisecorr_NVE.txt" using 1:($2/5.0) w lp title "T = 5.0" ls 1
replot "../OUTPUT/T2.0/N1080/noisecorr_NVE.txt" using 1:($2/2.0) w lp title "T = 2.0" ls 2
replot "../OUTPUT/T1.0/N1080/noisecorr_NVE.txt" using 1:($2/1.0) w lp title "T = 1.0" ls 3
replot "../OUTPUT/T0.6/N1080/noisecorr_NVE.txt" using 1:($2/0.6) w lp title "T = 0.6" ls 4
replot 0 ls 0 notitle


###############################################################
# Compare Noise Autocorrelation with Diagonal Autocorrelation #
###############################################################
reset

set title "Compare Unnormalized Noise Autocorrelation with Diagonal Autocorrelation (with factor 1/T)"
set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation"
load '~/.gnuplotting/gnuplot-palettes-master/paired.pal'
plot   "../OUTPUT/T5.0/N1080/noisecorr_NVE.txt" using 1:($2/5.0) w lp title "K(t) - T = 5.0" ls 1
replot "../OUTPUT/T5.0/N1080/Cd_NVE.txt" using ($1*0.0025):($2/5.0):($3/5.0) w errorl title "C_d(t) - T = 5.0" ls 2
replot "../OUTPUT/T2.0/N1080/noisecorr_NVE.txt" using 1:($2/2.0) w lp title "K(t) - T = 2.0" ls 3
replot "../OUTPUT/T2.0/N1080/Cd_NVE.txt" using ($1*0.0025):($2/2.0):($3/2.0) w errorl title "C_d(t) - T = 2.0" ls 4
replot "../OUTPUT/T1.0/N1080/noisecorr_NVE.txt" using 1:($2/1.0) w lp title "K(t) - T = 1.0" ls 5
replot "../OUTPUT/T1.0/N1080/Cd_NVE.txt" using ($1*0.0025):($2/1.0):($3/1.0) w errorl title "C_d(t) - T = 1.0" ls 6
replot "../OUTPUT/T0.6/N1080/noisecorr_NVE.txt" using 1:($2/0.6) w lp title "K(t) - T = 0.6" ls 7
replot "../OUTPUT/T0.6/N1080/Cd_NVE.txt" using ($1*0.0025):($2/0.6):($3/1.0) w errorl title "C_d(t) - T = 0.6" ls 8

set title "Compare Normalized Noise Autocorrelation with Diagonal Autocorrelation"
set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation"
load '~/.gnuplotting/gnuplot-palettes-master/paired.pal'
plot   "../OUTPUT/T5.0/N1080/noisecorr_NVE.txt" using 1:($3) w lp title "K(t) - T = 5.0" ls 1
norm5=system("awk  '(NR==2){print $2}'  ../OUTPUT/T5.0/N1080/Cd_NVE.txt")
replot "../OUTPUT/T5.0/N1080/Cd_NVE.txt" using ($1*0.0025):($2/norm5):($3/norm5) w errorl title "C_d(t) - T = 5.0" ls 2
replot "../OUTPUT/T2.0/N1080/noisecorr_NVE.txt" using 1:($3) w lp title "K(t) - T = 2.0" ls 3
norm2=system("awk  '(NR==2){print $2}'  ../OUTPUT/T2.0/N1080/Cd_NVE.txt")
replot "../OUTPUT/T2.0/N1080/Cd_NVE.txt" using ($1*0.0025):($2/norm2):($3/norm2) w errorl title "C_d(t) - T = 2.0" ls 4

replot "../OUTPUT/T1.0/N1080/noisecorr_NVE.txt" using 1:($3) w lp title "K(t) - T = 1.0" ls 5
norm1=system("awk  '(NR==2){print $2}'  ../OUTPUT/T1.0/N1080/Cd_NVE.txt")
replot "../OUTPUT/T1.0/N1080/Cd_NVE.txt" using ($1*0.0025):($2/norm1):($3/norm1) w errorl title "C_d(t) - T = 1.0" ls 6

replot "../OUTPUT/T0.6/N1080/noisecorr_NVE.txt" using 1:($3) w lp title "K(t) - T = 0.6" ls 7
norm06=system("awk  '(NR==2){print $2}'  ../OUTPUT/T0.6/N1080/Cd_NVE.txt")
replot "../OUTPUT/T0.6/N1080/Cd_NVE.txt" using ($1*0.0025):($2/norm06):($3/norm5) w errorl title "C_d(t) - T = 0.6" ls 8




########################################################################
# Compare different methods of evaluation of the Noise Autocorrelation #
########################################################################

set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation"
load '~/.gnuplotting/gnuplot-palettes-master/paired.pal'
plot   "../OUTPUT/T5.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "K^{NVE}" ls 1
replot   "../OUTPUT/T5.0/N1080/noisecorrlinear_NVE.txt" using 1:3 w lp title "K_{l}^{NVE}" ls 3
replot "../OUTPUT/T5.0/N1080/noisecorrlinear-selfcon_NVE.txt" using 1:3 w lp title "K_{lsc}^{NVE}" ls 5
replot 0 ls 0 notitle


set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation"
load '~/.gnuplotting/gnuplot-palettes-master/paired.pal'
plot   "../OUTPUT/T2.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "K^{NVE}" ls 1
replot   "../OUTPUT/T2.0/N1080/noisecorrlinear_NVE.txt" using 1:3 w lp title "K_{l}^{NVE}" ls 3
replot "../OUTPUT/T2.0/N1080/noisecorrlinear-selfcon_NVE.txt" using 1:3 w lp title "K_{lsc}^{NVE}" ls 5
replot 0 ls 0 notitle


set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation"
load '~/.gnuplotting/gnuplot-palettes-master/paired.pal'
plot   "../OUTPUT/T1.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "K^{NVE}" ls 1
replot   "../OUTPUT/T1.0/N1080/noisecorrlinear_NVE.txt" using 1:3 w lp title "K_{l}^{NVE}" ls 3
replot "../OUTPUT/T1.0/N1080/noisecorrlinear-selfcon_NVE.txt" using 1:3 w lp title "K_{lsc}^{NVE}" ls 5
replot 0 ls 0 notitle


