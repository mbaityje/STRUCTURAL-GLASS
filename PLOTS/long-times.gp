#!/bin/env gnuplot

# AT THE MOMENT I AM USING THE NORMALIZED CORRELATION FUNCTIONS, BUT I AM NOT SURE WHETHER I SHOULD NORMALIZE THEM


#################################################################
#                                                               #
# Fit the long-time behavior of the NOISE correlation functions #
#                                                               #
#################################################################

set logs x
set xlabel "{/Times-Italic t}"
set ylabel "{/Times-Italic K}({/Times-Italic t})"
plot[:][-.1:1.1]"../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 5.0" ls 1,\
"../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 2.0" ls 2,\
 "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 1.0" ls 3,\
 "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 0.8" ls 4,\
 "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 0.7" ls 5,\
 "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "T = 0.6" ls 6,\
 "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "T = 0.55" ls 7,\
 "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "T = 0.52" ls 8,\
 "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "T = 0.49" ls 9,\
 "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "T = 0.47" ls 10

f20(x)=A20*exp(-(x/t20)**b20)
f10(x)=A10*exp(-(x/t10)**b10)
f08(x)=A08*exp(-(x/t08)**b08)
f07(x)=A07*exp(-(x/t07)**b07)
f06(x)=A06*exp(-(x/t06)**b06)
f055(x)=A055*exp(-(x/t055)**b055)
f052(x)=A052*exp(-(x/t052)**b052)
f049(x)=A049*exp(-(x/t049)**b049)
f047(x)=A047*exp(-(x/t047)**b047); t047=10; b047=1; A047=.18
f046(x)=A046*exp(-(x/t046)**b046); t046=10; b046=1; A046=.18

set arrow nohead from 0.03,0 to 0.03,1
fit [0.4:] f046(x) "../OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via A046,t046,b046
plot[0.4:] "../OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.46" ls 11, f046(x)

fit [0.3:] f047(x) "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via A047,t047,b047
plot[0.3:] "../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.47" ls 10, f047(x)

fit [0.3:] f049(x) "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via A049,t049,b049
plot[0.3:] "../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.49" ls 10, f049(x)

fit [0.3:] f052(x) "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via A052,t052,b052
plot[0.3:] "../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.52" ls 10, f052(x)

fit [0.3:] f055(x) "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via A055,t055,b055
plot[0.3:] "../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.55" ls 10, f055(x)

fit [0.4:] f06(x) "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via A06,t06,b06
plot[0.4:] "../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.6" ls 10, f06(x)

fit [0.3:] f07(x) "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via A07,t07,b07
plot[0.3:] "../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.7" ls 10, f07(x)

fit [0.3:] f08(x) "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via A08,t08,b08
plot[0.3:] "../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 0.8" ls 10, f08(x)

fit [0.3:] f10(x) "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via A10,t10,b10
plot[0.3:] "../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 1.0" ls 10, f10(x)

fit [0.22:] f20(x) "../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M3.txt" u 1:3 via A20,t20,b20
plot[0.22:] "../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w errorl title "T = 2.0" ls 20, f20(x)


!rm -f "./SMALL-DATA/long-times_noise.txt"
set print "./SMALL-DATA/long-times_noise.txt"
pr "T plateau tau beta"
pr 2.0,A20,t20,b20
pr 1.0,A10,t10,b10
pr 0.8,A08,t08,b08
pr 0.7,A07,t07,b07
pr 0.6,A06,t06,b06
pr 0.55,A055,t055,b055
pr 0.52,A052,t052,b052
pr 0.49,A049,t049,b20
pr 0.47,A047,t047,b047
pr 0.46,A046,t046,b046

####################################################################
#                                                                  #
# Fit the long-time behavior of the DIAGONAL correlation functions #
#                                                                  #
####################################################################

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
plot[:][:]"../OUTPUT/T5.0/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm50):($3/norm50) w lp title "T = 5.0" ls 1,\
"../OUTPUT/T2.0/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm20):($3/norm20) w lp title "T = 2.0" ls 2,\
"../OUTPUT/T1.0/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm10):($3/norm10) w lp title "T = 1.0" ls 3,\
"../OUTPUT/T0.8/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm08):($3/norm08) w lp title "T = 0.8" ls 4,\
"../OUTPUT/T0.7/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm07):($3/norm07) w lp title "T = 0.7" ls 5,\
"../OUTPUT/T0.6/N1080/Cd_NVT.txt"  using ($1*0.0025):($2/norm06):($3/norm06) w lp title "T = 0.6" ls 6,\
"../OUTPUT/T0.55/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm055):($3/norm055) w lp title "T = 0.55" ls 7,\
"../OUTPUT/T0.52/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm052):($3/norm052) w lp title "T = 0.52" ls 8,\
"../OUTPUT/T0.49/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm049):($3/norm049) w lp title "T = 0.49" ls 9,\
"../OUTPUT/T0.47/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm047):($3/norm047) w lp title "T = 0.47" ls 10,\
"../OUTPUT/T0.46/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm046):($3/norm046) w lp title "T = 0.46" ls 11

set arrow nohead from 0.4,0 to 0.4,1
fit [0.4:] f046(x) "../OUTPUT/T0.46/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm50):($3/norm50) via A046,t046,b046
plot[:] "../OUTPUT/T0.46/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm50):($3/norm50) w errorl title "T = 0.46" ls 10, f046(x)

fit [0.4:] f047(x) "../OUTPUT/T0.47/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm50):($3/norm50) via A047,t047,b047
plot[:] "../OUTPUT/T0.47/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm50):($3/norm50) w errorl title "T = 0.47" ls 10, f047(x)

fit [0.4:] f049(x) "../OUTPUT/T0.49/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm50):($3/norm50) via A049,t049,b049
plot[:] "../OUTPUT/T0.49/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm50):($3/norm50) w errorl title "T = 0.49" ls 10, f049(x)

fit [0.4:] f052(x) "../OUTPUT/T0.52/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm50):($3/norm50) via A052,t052,b052
plot[:] "../OUTPUT/T0.52/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm50):($3/norm50) w errorl title "T = 0.52" ls 10, f052(x)

fit [0.4:] f055(x) "../OUTPUT/T0.55/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm50):($3/norm50) via A055,t055,b055
plot[:] "../OUTPUT/T0.55/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm50):($3/norm50) w errorl title "T = 0.55" ls 10, f055(x)

fit [0.4:] f06(x) "../OUTPUT/T0.6/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm50):($3/norm50) via A06,t06,b06
plot[:] "../OUTPUT/T0.6/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm50):($3/norm50) w errorl title "T = 0.6" ls 10, f06(x)

fit [0.4:] f07(x) "../OUTPUT/T0.7/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm50):($3/norm50) via A07,t07,b07
plot[:] "../OUTPUT/T0.7/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm50):($3/norm50) w errorl title "T = 0.7" ls 10, f07(x)

fit [0.4:] f08(x) "../OUTPUT/T0.8/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm50):($3/norm50) via A08,t08,b08
plot[:] "../OUTPUT/T0.8/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm50):($3/norm50) w errorl title "T = 0.8" ls 10, f08(x)

fit [0.4:] f10(x) "../OUTPUT/T1.0/N1080/Cd_NVT.txt" u ($1*0.0025):($2/norm50):($3/norm50) via A10,t10,b10
plot[:] "../OUTPUT/T1.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/norm50):($3/norm50) w errorl title "T = 1.0" ls 10, f10(x)


!rm -f "./SMALL-DATA/long-times_diag.txt"
set print "./SMALL-DATA/long-times_diag.txt"
pr "T plateau tau beta"
pr 2.0,A20,t20,b20
pr 1.0,A10,t10,b10
pr 0.8,A08,t08,b08
pr 0.7,A07,t07,b07
pr 0.6,A06,t06,b06
pr 0.55,A055,t055,b055
pr 0.52,A052,t052,b052
pr 0.49,A049,t049,b20
pr 0.47,A047,t047,b047
pr 0.46,A046,t046,b046



###################################
#                                 #
# Plot long-time behavior of both #
#                                 #
###################################

set term post enh c eps
set out "./FIGURES/plateau.eps"
reset
set ylabel "Plateau height"
set xlabel "{/Times-Italic T}"
p [0:2.5] "./SMALL-DATA/long-times_noise.txt" u 1:2:(0.06) w circles title "Noise",\
 "./SMALL-DATA/long-times_diag.txt" u 1:2:(0.06) w circles title "Diagonal"

set out "./FIGURES/beta.eps"
reset
set ylabel "Stretching exponent"
set xlabel "{/Times-Italic T}"
p [0:2.5] "./SMALL-DATA/long-times_noise.txt" u 1:4:(0.06) w circles title "Noise",\
"./SMALL-DATA/long-times_diag.txt" u 1:4:(0.06) w circles title "Diagonal"


reset
# Fit divergence of tau
Tc_noise=0.45; expo_noise=-2; amp_noise=1
Tc_diag =0.45; expo_diag =-2; amp_diag =1
fit [:1] amp_noise*(x-Tc_noise)**expo_noise "./SMALL-DATA/long-times_noise.txt" u 1:3 via amp_noise,expo_noise,Tc_noise
fit [:1] amp_diag*(x-Tc_diag)**expo_diag "./SMALL-DATA/long-times_diag.txt" u 1:3 via amp_diag,expo_diag,Tc_diag

set out "./FIGURES/tau.eps"
set ylabel "Relaxation time"
set xlabel "{/Times-Italic T}"
p [:2.5] "./SMALL-DATA/long-times_noise.txt" u 1:3 w p pt 11 ps 3 not, amp_noise*(x-Tc_noise)**expo_noise,\
"./SMALL-DATA/long-times_diag.txt" u 1:3 w p pt 11 ps 3 not, amp_noise*(x-Tc_noise)**expo_noise

set out "./FIGURES/tau-logs.eps"
reset
set logs
set xlabel "{/Times-Italic T-T}_c"
p "./SMALL-DATA/long-times_noise.txt" u ($1-Tc_noise):3 w p pt 11 ps 3, amp_noise*x**expo_noise,\
"./SMALL-DATA/long-times_diag.txt" u ($1-Tc_diag):3 w p pt 11 ps 3, amp_diag*x**expo_diag
