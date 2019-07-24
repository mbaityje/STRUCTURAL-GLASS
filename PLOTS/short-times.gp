#!/usr/bin/env gnuplot
reset
set term post enh c eps font "Times-Roman,32"
############################################
#                                          #
#  1/COSH FITS OF THE SHORT-TIME BEHAVIOR  #
#                                          #
############################################
top=1
bottom=3
left=6.1
right=1

set output "./FIGURES/short-times_cosh.eps"
set tmargin top
set bmargin bottom
set rmargin right
set lmargin left
set key bottom right
set xlabel "{/Times-Italic T}"
set ylabel "{/Times-Italic a}_2"
plot [0:5.1][0:]"< awk '(NR%3==0 && $1>=0.47)' ../OUTPUT/T*/N1080/shortNoise_NVT.txt"  using 1:5:6 w yerrorbars lc rgb "dark-violet" t "Noise",\
"< awk '(NR%3==0 && $1>=0.47)' ../OUTPUT/T*/N1080/shortDiag_NVT.txt"  using 1:5:6 w yerrorbars lc rgb "dark-green" t "Diagonal"
set out


#Plot correlation noise correlation function and its fit (I made a change in the columns but didnt check it's ok because I won't use these plots)
reset
set term post enh c eps font "Times-Roman,32"
set output "./FIGURES/short-times-fits_noise.eps"
set tmargin top
set bmargin bottom
set rmargin right
set lmargin left
set xtics 0,0.01
unset key
set xlabel "{/Times-Italic t}"
set ylabel "{/Times-Italic K}({/Times-Italic t})" offset 1.5,0
p[:0.04][0:] "../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 1,\
"../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 2,\
"../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 3,\
"../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 4,\
"../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 5,\
"../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 6,\
"../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 7,\
"../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 8,\
"../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 9,\
"../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 10,\
"../OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ls 11,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T5.0/N1080/shortNoise_NVT.txt`"*x) ls 1,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T2.0/N1080/shortNoise_NVT.txt`"*x) ls 2,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T1.0/N1080/shortNoise_NVT.txt`"*x) ls 3,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.8/N1080/shortNoise_NVT.txt`"*x) ls 4,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.7/N1080/shortNoise_NVT.txt`"*x) ls 5,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.6/N1080/shortNoise_NVT.txt`"*x) ls 6,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.55/N1080/shortNoise_NVT.txt`"*x) ls 7,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.52/N1080/shortNoise_NVT.txt`"*x) ls 8,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.49/N1080/shortNoise_NVT.txt`"*x) ls 9,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.47/N1080/shortNoise_NVT.txt`"*x) ls 10,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.46/N1080/shortNoise_NVT.txt`"*x) ls 11
set out

#Plot correlation diagonal correlation function and its fit
reset
set term post enh c eps font "Times-Roman,32"
set tmargin top
set bmargin bottom
set rmargin right
set lmargin left
set output "./FIGURES/short-times-fits_diag.eps"
unset key
set xtics 0,0.01
set xlabel "{/Times-Italic t}"
set ylabel "{/Times-Italic C}_d({/Times-Italic t})" offset 1.5,0
normCd50=`awk '($1==0){print $2}' ../OUTPUT/T5.0/N1080/Cd_NVT.txt`
normCd20=`awk '($1==0){print $2}' ../OUTPUT/T2.0/N1080/Cd_NVT.txt`
normCd10=`awk '($1==0){print $2}' ../OUTPUT/T1.0/N1080/Cd_NVT.txt`
normCd08=`awk '($1==0){print $2}' ../OUTPUT/T0.8/N1080/Cd_NVT.txt`
normCd07=`awk '($1==0){print $2}' ../OUTPUT/T0.7/N1080/Cd_NVT.txt`
normCd06=`awk '($1==0){print $2}' ../OUTPUT/T0.6/N1080/Cd_NVT.txt`
normCd055=`awk '($1==0){print $2}' ../OUTPUT/T0.55/N1080/Cd_NVT.txt`
normCd052=`awk '($1==0){print $2}' ../OUTPUT/T0.52/N1080/Cd_NVT.txt`
normCd049=`awk '($1==0){print $2}' ../OUTPUT/T0.49/N1080/Cd_NVT.txt`
normCd047=`awk '($1==0){print $2}' ../OUTPUT/T0.47/N1080/Cd_NVT.txt`
normCd046=`awk '($1==0){print $2}' ../OUTPUT/T0.46/N1080/Cd_NVT.txt`
p[:0.04][0:1.] "../OUTPUT/T5.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd50):($3/normCd50) w errorl notitle ls 1 pt 2,\
"../OUTPUT/T2.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd20):($3/normCd20) w errorl notitle ls 2 pt 2,\
"../OUTPUT/T1.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd10):($3/normCd10) w errorl notitle ls 3 pt 2,\
"../OUTPUT/T0.8/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd08):($3/normCd08) w errorl notitle ls 4 pt 2,\
"../OUTPUT/T0.7/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd07):($3/normCd07) w errorl notitle ls 5 pt 2,\
"../OUTPUT/T0.6/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd06):($3/normCd06) w errorl notitle ls 6 pt 2,\
"../OUTPUT/T0.55/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd055):($3/normCd055) w errorl notitle ls 7 pt 2,\
"../OUTPUT/T0.52/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd052):($3/normCd052) w errorl notitle ls 8 pt 2,\
"../OUTPUT/T0.49/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd049):($3/normCd049) w errorl notitle ls 9 pt 2,\
"../OUTPUT/T0.47/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd047):($3/normCd047) w errorl notitle ls 10 pt 2,\
"../OUTPUT/T0.46/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd046):($3/normCd046) w errorl notitle ls 11 pt 2,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T5.0/N1080/shortDiag_NVT.txt`"*x) ls 1,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T2.0/N1080/shortDiag_NVT.txt`"*x) ls 2,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T1.0/N1080/shortDiag_NVT.txt`"*x) ls 3,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.8/N1080/shortDiag_NVT.txt`"*x) ls 4,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.7/N1080/shortDiag_NVT.txt`"*x) ls 5,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.6/N1080/shortDiag_NVT.txt`"*x) ls 6,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.55/N1080/shortDiag_NVT.txt`"*x) ls 7,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.52/N1080/shortDiag_NVT.txt`"*x) ls 8,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.49/N1080/shortDiag_NVT.txt`"*x) ls 9,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.47/N1080/shortDiag_NVT.txt`"*x) ls 10,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.46/N1080/shortDiag_NVT.txt`"*x) ls 11
set out



set output "./FIGURES/short-times_inset-collapse.eps"
set multiplot
set tmargin top
set bmargin bottom
set rmargin right
set lmargin left
set ytics 0,0.2
set xtics format "10^{%T}"
set xlabel "{/Times-Italic t}"
set ylabel "Autocorrelation"
set logs x
p[0.001:1][0:] 1./cosh("`awk '(NR==3){print $4}' ../OUTPUT/T5.0/N1080/shortNoise_NVT.txt`"*x) ls 1,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T2.0/N1080/shortNoise_NVT.txt`"*x) ls 2,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T1.0/N1080/shortNoise_NVT.txt`"*x) ls 3,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.8/N1080/shortNoise_NVT.txt`"*x) ls 4,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.7/N1080/shortNoise_NVT.txt`"*x) ls 5,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.6/N1080/shortNoise_NVT.txt`"*x) ls 6,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.55/N1080/shortNoise_NVT.txt`"*x) ls 7,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.52/N1080/shortNoise_NVT.txt`"*x) ls 8,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.49/N1080/shortNoise_NVT.txt`"*x) ls 9,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.47/N1080/shortNoise_NVT.txt`"*x) ls 10,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.46/N1080/shortNoise_NVT.txt`"*x) ls 11,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T5.0/N1080/shortDiag_NVT.txt`"*x) ls 1,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T2.0/N1080/shortDiag_NVT.txt`"*x) ls 2,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T1.0/N1080/shortDiag_NVT.txt`"*x) ls 3,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.8/N1080/shortDiag_NVT.txt`"*x) ls 4,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.7/N1080/shortDiag_NVT.txt`"*x) ls 5,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.6/N1080/shortDiag_NVT.txt`"*x) ls 6,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.55/N1080/shortDiag_NVT.txt`"*x) ls 7,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.52/N1080/shortDiag_NVT.txt`"*x) ls 8,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.49/N1080/shortDiag_NVT.txt`"*x) ls 9,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.47/N1080/shortDiag_NVT.txt`"*x) ls 10,\
1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.46/N1080/shortDiag_NVT.txt`"*x) ls 11
set tmargin 2
set bmargin 10
set rmargin 2
set lmargin 22
set ytics 0,0.3 font ",24"
set xtics format "10^{%T}" font ",24"
unset ylabel
unset xlabel
rescale=1.15
#p[0.001:1][0:] 1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T5.0/N1080/shortNoise_NVT.txt`"*x) ls 1,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T2.0/N1080/shortNoise_NVT.txt`"*x) ls 2,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T1.0/N1080/shortNoise_NVT.txt`"*x) ls 3,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.8/N1080/shortNoise_NVT.txt`"*x) ls 4,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.7/N1080/shortNoise_NVT.txt`"*x) ls 5,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.6/N1080/shortNoise_NVT.txt`"*x) ls 6,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.55/N1080/shortNoise_NVT.txt`"*x) ls 7,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.52/N1080/shortNoise_NVT.txt`"*x) ls 8,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.49/N1080/shortNoise_NVT.txt`"*x) ls 9,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.47/N1080/shortNoise_NVT.txt`"*x) ls 10,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.46/N1080/shortNoise_NVT.txt`"*x) ls 11,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T5.0/N1080/shortDiag_NVT.txt`"*x*rescale) ls 1,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T2.0/N1080/shortDiag_NVT.txt`"*x*rescale) ls 2,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T1.0/N1080/shortDiag_NVT.txt`"*x*rescale) ls 3,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.8/N1080/shortDiag_NVT.txt`"*x*rescale) ls 4,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.7/N1080/shortDiag_NVT.txt`"*x*rescale) ls 5,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.6/N1080/shortDiag_NVT.txt`"*x*rescale) ls 6,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.55/N1080/shortDiag_NVT.txt`"*x*rescale) ls 7,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.52/N1080/shortDiag_NVT.txt`"*x*rescale) ls 8,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.49/N1080/shortDiag_NVT.txt`"*x*rescale) ls 9,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.47/N1080/shortDiag_NVT.txt`"*x*rescale) ls 10,\
#1./cosh("`awk '(NR==3){print $5}' ../OUTPUT/T0.46/N1080/shortDiag_NVT.txt`"*x*rescale) ls 11
unset multiplot
set out



############################################
#                                          #
# GAUSSIAN FITS OF THE SHORT-TIME BEHAVIOR #
#                                          #
############################################
set term unknown
set xlabel "t"
set ylabel "K(t)"
plot[:0.05][-.1:1.1]"../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 5.0" ls 1
replot "../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 2.0" ls 2
replot "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 1.0" ls 3
replot "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 0.8" ls 4
replot "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 0.7" ls 5
replot "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 0.6" ls 6
replot "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp title "T = 0.55" ls 7
replot "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp title "T = 0.52" ls 8
replot "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp title "T = 0.49" ls 9
replot "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp title "T = 0.47" ls 10
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
fit [:0.03] f047(x) "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via a047
plot "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.47" ls 10, f047(x)

fit [:0.03] f049(x) "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via a049
plot "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.49" ls 10, f049(x)

fit [:0.03] f052(x) "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via a052
plot "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.52" ls 10, f052(x)

fit [:0.03] f055(x) "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via a055
plot "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.55" ls 10, f055(x)

fit [:0.03] f06(x) "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via a06
plot "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.6" ls 10, f06(x)

fit [:0.03] f07(x) "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via a07
plot "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.7" ls 10, f07(x)

fit [:0.03] f08(x) "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via a08
plot "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.8" ls 10, f08(x)

fit [:0.03] f10(x) "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via a10
plot "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 1.0" ls 10, f10(x)

fit [:0.015] f20(x) "../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via a20
plot "../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 2.0" ls 20, f20(x)

fit [:0.01] f50(x) "../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via a50
plot "../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 5.0" ls 50, f50(x)

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
