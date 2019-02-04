#!/usr/bin/env gnuplot
reset
set term post enh c eps #font "Times-Roman,32"
############################################
#                                          #
#  1/COSH FITS OF THE SHORT-TIME BEHAVIOR  #
#                                          #
############################################
set output "./FIGURES/short-times_cosh.eps"
set xlabel "{/Times-Italic t}"
set ylabel "{/Times-Italic a}_2"
plot [0:][0:]"< awk '(NR%3==0 && $1>=0.47)' ../OUTPUT/T*/N1080/shortNoise_NVT.txt"  using 1:4:(0.06) w circles t "Noise",\
"< awk '(NR%3==0 && $1>=0.47)' ../OUTPUT/T*/N1080/shortDiag_NVT.txt"  using 1:4:(0.06) w circles t "Diagonal"

plot [0:][0:]"< awk '(NR%3==0 && $1>=0.47)' ../OUTPUT/T*/N1080/shortNoise_NVT.txt"  using 1:($4):(0.06) w circles t "Noise",\
"< awk '(NR%3==0 && $1>=0.47)' ../OUTPUT/T*/N1080/shortDiag_NVT.txt"  using 1:($4*1.15):(0.06) w circles t "Diagonal x 1.15"
set out





############################################
#                                          #
# GAUSSIAN FITS OF THE SHORT-TIME BEHAVIOR #
#                                          #
############################################
set term unknown
set xlabel "t"
set ylabel "K(t)"
plot[:0.05][-.1:1.1]"../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 5.0" ls 1
replot "../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 2.0" ls 2
replot "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 1.0" ls 3
replot "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 0.8" ls 4
replot "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 0.7" ls 5
replot "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 0.6" ls 6
replot "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "T = 0.55" ls 7
replot "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "T = 0.52" ls 8
replot "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "T = 0.49" ls 9
replot "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "T = 0.47" ls 10
replot 0 ls 0 notitle

f50(x)=1-a50*x**2
f20(x)=1-a20*x**2
f10(x)=1-a10*x**2
f08(x)=1-a08*x**2
f07(x)=1-a07*x**2
f06(x)=1-a06*x**2
f055(x)=1-a055*x**2
f052(x)=1-a052*x**2
f049(x)=1-a049*x**2
f047(x)=1-a047*x**2

set xr[0:0.04]
set arrow nohead from 0.03,0 to 0.03,1
fit [:0.03] f047(x) "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via a047
plot "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.47" ls 10, f047(x)

fit [:0.03] f049(x) "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via a049
plot "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.49" ls 10, f049(x)

fit [:0.03] f052(x) "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via a052
plot "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.52" ls 10, f052(x)

fit [:0.03] f055(x) "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via a055
plot "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.55" ls 10, f055(x)

fit [:0.03] f06(x) "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via a06
plot "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.6" ls 10, f06(x)

fit [:0.03] f07(x) "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via a07
plot "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.7" ls 10, f07(x)

fit [:0.03] f08(x) "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via a08
plot "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.8" ls 10, f08(x)

fit [:0.03] f10(x) "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via a10
plot "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 1.0" ls 10, f10(x)

fit [:0.015] f20(x) "../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via a20
plot "../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 2.0" ls 20, f20(x)

fit [:0.01] f50(x) "../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via a50
plot "../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 5.0" ls 50, f50(x)

!rm -f "./SMALL-DATA/short-times_gauss.txt"
set print "./SMALL-DATA/short-times_gauss.txt"
pr "T alpha"
pr 5.0,a50
pr 2.0,a20
pr 1.0,a10
pr 0.8,a08
pr 0.7,a07
pr 0.6,a06
pr 0.55,a055
pr 0.52,a052
pr 0.49,a049
pr 0.47,a047



set term post enh c eps
set out "./FIGURES/short-times_gauss.eps"
reset
set ylabel "{/Symbol a}"
set xlabel "{/Times-Italic T}"
p [0:2.5] "./SMALL-DATA/short-times_gauss.txt" u 1:2:(0.06) w circles notitle
