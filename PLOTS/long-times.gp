#!/bin/env gnuplot

set fit errorvariables

#################################################################
#                                                               #
# Fit the long-time behavior of the NOISE correlation functions #
#                                                               #
#################################################################

set logs x
set xlabel "{/Times-Italic t}"
set ylabel "{/Times-Italic K}({/Times-Italic t})"
plot[:][-.1:1.1]"../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 5.0" ls 1,\
"../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 2.0" ls 2,\
 "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 1.0" ls 3,\
 "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 0.8" ls 4,\
 "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 0.7" ls 5,\
 "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp title "T = 0.6" ls 6,\
 "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp title "T = 0.55" ls 7,\
 "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp title "T = 0.52" ls 8,\
 "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp title "T = 0.49" ls 9,\
 "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp title "T = 0.47" ls 10

f50(x)=A50*exp(-(x/t50)**b50)
f20(x)=A20*exp(-(x/t20)**b20)
f10(x)=A10*exp(-(x/t10)**b10)
f08(x)=A08*exp(-(x/t08)**b08)
f07(x)=A07*exp(-(x/t07)**b07); t07=10; b07=1; A07=.18
f06(x)=A06*exp(-(x/t06)**b06); t06=10; b06=1; A06=.18
f055(x)=A055*exp(-(x/t055)**b055); t055=10; b055=1; A055=.18
f052(x)=A052*exp(-(x/t052)**b052); t052=10; b052=1; A052=.18
f049(x)=A049*exp(-(x/t049)**b049); t049=10; b049=1; A049=.18
f047(x)=A047*exp(-(x/t047)**b047); t047=10; b047=1; A047=.18
f046(x)=A046*exp(-(x/t046)**b046); t046=10; b046=1; A046=.18
f045(x)=A045*exp(-(x/t045)**b045); t045=10; b046=1; A045=.18

unset arrow
left=20
set arrow nohead from left,0 to left,10
fit [20:] f045(x) "../OUTPUT/T0.45/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A045,t045,b045
plot[:] "../OUTPUT/T0.45/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.45" ls 12, f045(x)

fit [1:] f046(x) "../OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A046,t046,b046
plot[:] "../OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.46" ls 11, f046(x)

fit [1:] f047(x) "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A047,t047,b047
plot[:] "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.47" ls 10, f047(x)

fit [1:] f049(x) "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A049,t049,b049
plot[:] "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.49" ls 9, f049(x)

fit [1:] f052(x) "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A052,t052,b052
plot[:] "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.52" ls 8, f052(x)

fit [1:] f055(x) "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A055,t055,b055
plot[:] "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.55" ls 7, f055(x)

fit [1:] f06(x) "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A06,t06,b06
plot[:] "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.6" ls 6, f06(x)

fit [1:] f07(x) "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A07,t07,b07
plot[:] "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.7" ls 5, f07(x)

fit [1:] f08(x) "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A08,t08,b08
plot[:] "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 0.8" ls 4, f08(x)

fit [.6:] f10(x) "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A10,t10,b10
plot[:] "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 1.0" ls 3, f10(x)

fit [.3:] f20(x) "../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 via A20,t20,b20
plot[:] "../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w errorl title "T = 2.0" ls 2, f20(x)


!rm -f "./SMALL-DATA/long-times_noise.txt"
set print "./SMALL-DATA/long-times_noise.txt"
pr "T plateau perr tau terr beta berr"
pr 2.0,   A20, A20_err,  t20, t20_err,  b20,  b20_err
pr 1.0,   A10, A10_err,  t10, t10_err,  b10,  b10_err
pr 0.8,   A08, A08_err,  t08, t08_err,  b08,  b08_err
pr 0.7,   A07, A07_err,  t07, t07_err,  b07,  b07_err
pr 0.6,   A06, A06_err,  t06, t06_err,  b06,  b06_err
pr 0.55, A055,A055_err, t055,t055_err, b055, b055_err
pr 0.52, A052,A052_err, t052,t052_err, b052, b052_err
pr 0.49, A049,A049_err, t049,t049_err, b049, b049_err
pr 0.47, A047,A047_err, t047,t047_err, b047, b047_err
pr 0.46, A046,A046_err, t046,t046_err, b046, b046_err
pr 0.45, A045,A045_err, t045,t045_err, b045, b045_err
set print

####################################################################
#                                                                  #
# Fit the long-time behavior of the DIAGONAL correlation functions #
#                                                                  #
####################################################################
unset arrow
set xr[*:*]; set yr[*:*]
set logs x
set xlabel "{/Times-Italic t}"
set ylabel "{/Times-Italic C}_d({/Times-Italic t})"
norm50="`awk '($1==0){print $2}' ../OUTPUT/T5.0/N1080/Cd_NVT.txt`"
norm20="`awk '($1==0){print $2}' ../OUTPUT/T2.0/N1080/Cd_NVT.txt`"
norm10="`awk '($1==0){print $2}' ../OUTPUT/T1.0/N1080/Cd_NVT.txt`"
norm08="`awk '($1==0){print $2}' ../OUTPUT/T0.8/N1080/Cd_NVT.txt`"
norm07="`awk '($1==0){print $2}' ../OUTPUT/T0.7/N1080/Cd_NVT.txt`"
norm06="`awk '($1==0){print $2}' ../OUTPUT/T0.6/N1080/Cd_NVT.txt`"
norm055="`awk '($1==0){print $2}' ../OUTPUT/T0.55/N1080/Cd_NVT.txt`"
norm052="`awk '($1==0){print $2}' ../OUTPUT/T0.52/N1080/Cd_NVT.txt`"
norm049="`awk '($1==0){print $2}' ../OUTPUT/T0.49/N1080/Cd_NVT.txt`"
norm047="`awk '($1==0){print $2}' ../OUTPUT/T0.47/N1080/Cd_NVT.txt`"
norm046="`awk '($1==0){print $2}' ../OUTPUT/T0.46/N1080/Cd_NVT.txt`"
norm045="`awk '($1==0){print $2}' ../OUTPUT/T0.45/N1080/Cd_NVT.txt`"
plot[:][:]"../OUTPUT/T5.0/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm50):($3/norm50) w lp ls 1 title "T = 5.0",\
"../OUTPUT/T2.0/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm20):($3/norm20) w lp title "T = 2.0" ls 2,\
"../OUTPUT/T1.0/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm10):($3/norm10) w lp title "T = 1.0" ls 3,\
"../OUTPUT/T0.8/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm08):($3/norm08) w lp title "T = 0.8" ls 4,\
"../OUTPUT/T0.7/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm07):($3/norm07) w lp title "T = 0.7" ls 5,\
"../OUTPUT/T0.6/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm06):($3/norm06) w lp title "T = 0.6" ls 6,\
"../OUTPUT/T0.55/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm055):($3/norm055) w lp title "T = 0.55" ls 7,\
"../OUTPUT/T0.52/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm052):($3/norm052) w lp title "T = 0.52" ls 8,\
"../OUTPUT/T0.49/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm049):($3/norm049) w lp title "T = 0.49" ls 9,\
"../OUTPUT/T0.47/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm047):($3/norm047) w lp title "T = 0.47" ls 10,\
"../OUTPUT/T0.46/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm046):($3/norm046) w lp title "T = 0.46" ls 11,\
"../OUTPUT/T0.45/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm045):($3/norm045) w lp title "T = 0.45" ls 12

set arrow nohead from 1,0 to 1,1
fit [3:] f045(x) "../OUTPUT/T0.45/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm045):($3/norm045) yerror via A045,t045,b045
plot[:] "../OUTPUT/T0.45/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm045):($3/norm045) w errorl title "T = 0.45" ls 10, f045(x)

fit [1:] f046(x) "../OUTPUT/T0.46/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm046):($3/norm046) yerror via A046,t046,b046
plot[:] "../OUTPUT/T0.46/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm046):($3/norm046) w errorl title "T = 0.46" ls 10, f046(x)

fit [1:] f047(x) "../OUTPUT/T0.47/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm047):($3/norm047) yerror via A047,t047,b047
plot[:] "../OUTPUT/T0.47/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm047):($3/norm047) w errorl title "T = 0.47" ls 10, f047(x)

fit [1:] f049(x) "../OUTPUT/T0.49/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm049):($3/norm049) yerror via A049,t049,b049
plot[:] "../OUTPUT/T0.49/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm049):($3/norm049) w errorl title "T = 0.49" ls 10, f049(x)

fit [1:] f052(x) "../OUTPUT/T0.52/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm052):($3/norm052) yerror via A052,t052,b052
plot[:] "../OUTPUT/T0.52/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm052):($3/norm052) w errorl title "T = 0.52" ls 10, f052(x)

fit [0.4:] f055(x) "../OUTPUT/T0.55/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm055):($3/norm055) yerror via A055,t055,b055
plot[:] "../OUTPUT/T0.55/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm055):($3/norm055) w errorl title "T = 0.55" ls 10, f055(x)

fit [0.4:] f06(x) "../OUTPUT/T0.6/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm06):($3/norm06) yerror via A06,t06,b06
plot[:] "../OUTPUT/T0.6/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm06):($3/norm06) w errorl title "T = 0.6" ls 10, f06(x)

fit [0.4:] f07(x) "../OUTPUT/T0.7/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm07):($3/norm07) yerror via A07,t07,b07
plot[:] "../OUTPUT/T0.7/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm07):($3/norm07) w errorl title "T = 0.7" ls 10, f07(x)

fit [0.4:] f08(x) "../OUTPUT/T0.8/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm08):($3/norm08) yerror via A08,t08,b08
plot[:] "../OUTPUT/T0.8/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm08):($3/norm08) w errorl title "T = 0.8" ls 10, f08(x)

fit [0.4:] f10(x) "../OUTPUT/T1.0/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm10):($3/norm10) yerror via A10,t10,b10
plot[:] "../OUTPUT/T1.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm10):($3/norm10) w errorl title "T = 1.0" ls 10, f10(x)

fit [0.4:] f20(x) "../OUTPUT/T2.0/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm20):($3/norm20) yerror via A20,t20,b20
plot[:] "../OUTPUT/T2.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm20):($3/norm20) w errorl title "T = 2.0" ls 10, f20(x)

fit [0.4:] f50(x) "../OUTPUT/T5.0/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm50):($3/norm50) yerror via A50,t50,b50
plot[:] "../OUTPUT/T5.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm50):($3/norm50) w errorl title "T = 5.0" ls 10, f50(x)

!rm -f "./SMALL-DATA/long-times_diag.txt"
set print "./SMALL-DATA/long-times_diag.txt"
pr "T plateau perr tau terr beta berr"
pr  5.0,  A50, A50_err,  t50, t50_err,  b50, b50_err
pr  2.0,  A20, A20_err,  t20, t20_err,  b20, b20_err
pr  1.0,  A10, A10_err,  t10, t10_err,  b10, b10_err
pr  0.8,  A08, A08_err,  t08, t08_err,  b08, b08_err
pr  0.7,  A07, A07_err,  t07, t07_err,  b07, b07_err
pr  0.6,  A06, A06_err,  t06, t06_err,  b06, b06_err
pr 0.55, A055,A055_err, t055,t055_err, b055,b055_err
pr 0.52, A052,A052_err, t052,t052_err, b052,b052_err
pr 0.49, A049,A049_err, t049,t049_err, b049,b049_err
pr 0.47, A047,A047_err, t047,t047_err, b047,b047_err
pr 0.46, A046,A046_err, t046,t046_err, b046,b046_err
pr 0.45, A045,A045_err, t045,t045_err, b045,b045_err
set print


###################################
#                                 #
# Plot long-time behavior of both #
#                                 #
###################################
reset
left=6
right=1
top=1

set term post enh c eps font "Times-Roman,32"
set out "./FIGURES/plateau.eps"
set lmargin left
set rmargin right
set tmargin top
set ylabel "Plateau height" offset 2,0
set xlabel "{/Times-Italic T}" offset 0,1
set ytics 0,0.1 font ",28" offset 0.5,0
set xtics 0,0.5 font ",28" offset 0,0.5
set key center right invert
p [0:2.1][0:] "./SMALL-DATA/long-times_noise.txt" u 1:2:3 w errorl title "Noise",\
 "./SMALL-DATA/long-times_diag.txt" u 1:2:3 w errorl title "Diagonal"

set out "./FIGURES/beta.eps"
set lmargin left
set rmargin right
set tmargin top
set ylabel "Stretching exponent" offset 2,0
set xlabel "{/Times-Italic T}" offset 0,1
set ytics 0,0.1 font ",28" offset 0.5,0
set xtics 0,0.5 font ",28" offset 0,0.5
set key center right noinvert
p [0:2.1][.5:1] "./SMALL-DATA/long-times_noise.txt" u 1:6:7 w errorl title "Noise",\
"./SMALL-DATA/long-times_diag.txt" u 1:6:7 w errorl title "Diagonal"


# Fit divergence of tau
Tc_noise=0.435; expo_noise=-2; amp_noise=1
Tc_diag =0.435; expo_diag =-2; amp_diag =1
Tc_MCT  =0.435; expo_MCT  =-2; amp_MCT  =1
maxfit=0.7
fit [:maxfit] amp_noise*(x-Tc_noise)**expo_noise "./SMALL-DATA/long-times_noise.txt" u 1:4:5 via amp_noise,expo_noise,Tc_noise
fit [:maxfit] amp_diag*(x-Tc_diag)**expo_diag "./SMALL-DATA/long-times_diag.txt" u 1:4:5 via amp_diag,expo_diag,Tc_diag
fit [:maxfit] amp_MCT*(x-Tc_MCT)**expo_MCT "../THERMALIZE/data/tauMCT.txt" u 1:2:3 via amp_MCT,expo_MCT,Tc_MCT

set out "./FIGURES/tau.eps"
set lmargin left
set rmargin right
set tmargin top
set ylabel "Relaxation time" offset 2,0
set xlabel "{/Times-Italic T}" offset 0,1
set ytics auto font ",28" offset 0.5,0
set xtics 0,0.5 font ",28" offset 0,0.5
p [:1] "./SMALL-DATA/long-times_noise.txt" u 1:4:5 w p pt 11 ps 3 not, amp_noise*(x-Tc_noise)**expo_noise not,\
"./SMALL-DATA/long-times_diag.txt" u 1:4:5 w p pt 11 ps 3 not, amp_diag*(x-Tc_diag)**expo_diag not

set out "./FIGURES/tau-logs.eps"
set lmargin left
set rmargin right
set tmargin top
reset
set logs
set xlabel "{/Times-Italic T-T}_{fit}"
p "./SMALL-DATA/long-times_noise.txt" u ($1-Tc_noise):4:5 w errorl pt 11 ps 1, amp_noise*x**expo_noise,\
"./SMALL-DATA/long-times_diag.txt" u ($1-Tc_diag):4:5 w errorl pt 11 ps 1, amp_diag*x**expo_diag




reset
set out "./FIGURES/tau_inset.eps"
set multiplot
set lmargin at screen 0.18
set rmargin at screen 0.98
set tmargin at screen 0.975
set bmargin at screen 0.13
set ylabel "{/Symbol-Oblique t}" offset 2,0
set xlabel "{/Times-Italic T}" offset 0,1
set ytics auto font ",28" offset 0.5,0 nomirror
set xtics auto font ",28" offset 0,0.5 nomirror
unset key
p [:1] "./SMALL-DATA/long-times_noise.txt" u 1:4:5 w errorl lc rgb "dark-violet" pt 7 ps 2 not,\
"./SMALL-DATA/long-times_diag.txt" u 1:4:5 w errorl lc rgb "dark-green" pt 6 ps 2 not,\
"../THERMALIZE/data/tauMCT.txt" u ($1):2:3 w errorl lc rgb "dark-orange" pt 4 ps 2 not
set lmargin at screen 0.45
set rmargin at screen 0.98
set tmargin at screen 0.975
set bmargin at screen 0.35
set logs
set xlabel "{/Times-Italic T-T}_{fit}"
unset ylabel
set ytics auto font ",24" offset 0.5,0 format "10^{%T}"
set xtics auto font ",24" offset 0,0.5
set key top right samplen 1
p[][:4000] "./SMALL-DATA/long-times_diag.txt" u ($1-Tc_diag):4:5 w errorl lc rgb "dark-green" pt 6 ps 1.3 t"{/Times-Italic C}_d({/Times-Italic t})", (x<maxfit-Tc_diag ? amp_diag*x**expo_diag : 1/0) not lc -1 lw .1,\
"../THERMALIZE/data/tauMCT.txt" u ($1-Tc_MCT):2:3 w errorl lc rgb "dark-orange" pt 4 ps 1.3  t"{/Times-Italic K}_{MCT}({/Times-Italic t})", (x<maxfit-Tc_MCT ? amp_MCT*x**expo_MCT : 1/0) not lc -1 lw .1,\
"./SMALL-DATA/long-times_noise.txt" u ($1-Tc_noise):4:5 w errorl lc rgb "dark-violet" pt 7 ps 1.3 t"{/Times-Italic K}({/Times-Italic t})", (x<maxfit-Tc_noise ? amp_noise*x**expo_noise : 1/0) not lc -1 lw .1
unset multiplot

pr "Tc_noise = ",Tc_noise," pm ", Tc_noise_err
pr "Tc_diag = ",Tc_diag," pm ", Tc_diag_err
pr "Tc_MCT = ",Tc_MCT," pm ", Tc_MCT_err
pr "expo_noise = ",expo_noise," pm ", expo_noise_err
pr "expo_diag = ",expo_diag," pm ", expo_diag_err
pr "expo_MCT = ",expo_MCT," pm ", expo_MCT_err
