#!/usr/bin/env gnuplot

#kmax for Kmct:
kmax="_kmax28"
#kmax=""


reset
set term post enh c eps font "Times-Roman,32"

########################################
# Show noise autocorrelation functions #
########################################

set output "./FIGURES/K-Cd.eps"
set lmargin 6
set tmargin .5
set rmargin .5
set bmargin 2.5
set logs x
set xlabel "{/Times-Italic t}" offset 0,.8
set ylabel "Memory ({/Times-Italic t})" offset 2.6,0
set key invert font ",27" samplen 1 maxrows 6
set tics font ",24"
set xtics format "10^{%T}" offset 0,0.4
set ytics offset 0.5,0
set label "{/Times-Italic K}({/Times-Italic t})" at 0.005,0.1
set label "{/Times-Italic C}_d({/Times-Italic t})" at 100,0.35
load 'scl.pal'
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
normCd045=`awk '($1==0){print $2}' ../OUTPUT/T0.45/N1080/Cd_NVT.txt`
plot[][-.1:1.1]"../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 5.0  " ls 1 pt 1 ps .8,\
"../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 2.0  " ls 2 pt 1 ps .8,\
"../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 1.0  " ls 3 pt 1 ps .8,\
"../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 0.8  " ls 4 pt 1 ps .8,\
"../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 0.7  " ls 5 pt 1 ps .8,\
"../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 0.6  " ls 6 pt 1 ps .8,\
"../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.55" ls 7 pt 1 ps .8,\
"../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.52" ls 8 pt 1 ps .8,\
"../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.49" ls 9 pt 1 ps .8,\
"../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.47" ls 10 pt 1 ps .8,\
"../OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.46" ls 11 pt 1 ps .8,\
"../OUTPUT/T0.45/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.45" ls 12 pt 1 ps .8,\
0 ls 0 notitle,\
"../OUTPUT/T5.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd50):($3/normCd50) w errorl  title "{/Times-Italic T} = 5.0  " ls 1 pt 2,\
"../OUTPUT/T2.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd20):($3/normCd20) w errorl title "{/Times-Italic T} = 2.0  " ls 2 pt 2,\
"../OUTPUT/T1.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd10):($3/normCd10) w errorl title "{/Times-Italic T} = 1.0  " ls 3 pt 2,\
"../OUTPUT/T0.8/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd08):($3/normCd08) w errorl title "{/Times-Italic T} = 0.8  " ls 4 pt 2,\
"../OUTPUT/T0.7/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd07):($3/normCd07) w errorl title "{/Times-Italic T} = 0.7  " ls 5 pt 2,\
"../OUTPUT/T0.6/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd06):($3/normCd06) w errorl title "{/Times-Italic T} = 0.6  " ls 6 pt 2,\
"../OUTPUT/T0.55/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd055):($3/normCd055) w errorl title "{/Times-Italic T} = 0.55" ls 7 pt 2,\
"../OUTPUT/T0.52/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd052):($3/normCd052) w errorl title "{/Times-Italic T} = 0.52" ls 8 pt 2,\
"../OUTPUT/T0.49/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd049):($3/normCd049) w errorl title "{/Times-Italic T} = 0.49" ls 9 pt 2,\
"../OUTPUT/T0.47/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd047):($3/normCd047) w errorl title "{/Times-Italic T} = 0.47" ls 10 pt 2,\
"../OUTPUT/T0.46/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd046):($3/normCd046) w errorl title "{/Times-Italic T} = 0.46" ls 11 pt 2,\
"../OUTPUT/T0.45/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd045):($3/normCd045) w errorl title "{/Times-Italic T} = 0.45" ls 12 pt 2
set out
unset label


# Only noise correlation function, without normalizing, and including factor 1/T
reset
set output "./FIGURES/KonT_unnorm.eps"
set xlabel "time [LJ units]"
set ylabel "K(t)/T" offset 2.5,0
set key invert font ",27" samplen 2 maxrows 5
set tics font ",24"
set xtics format "10^{%T}"
load 'scl.pal'
set logs x
plot   "../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/5.0) w lp title "{/Times-Italic T} = 5.0  " ls 1,\
"../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/2.0) w lp title "{/Times-Italic T} = 2.0  " ls 2,\
"../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/1.0) w lp title "{/Times-Italic T} = 1.0  " ls 3,\
"../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/0.8) w lp title "{/Times-Italic T} = 0.8  " ls 4,\
"../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/0.7) w lp title "{/Times-Italic T} = 0.7  " ls 5,\
"../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/0.6) w lp title "{/Times-Italic T} = 0.6  " ls 6,\
"../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/0.55) w lp title "{/Times-Italic T} = 0.55" ls 7,\
"../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/0.52) w lp title "{/Times-Italic T} = 0.52" ls 8,\
"../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/0.49) w lp title "{/Times-Italic T} = 0.49" ls 9,\
"../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/0.47) w lp title "{/Times-Italic T} = 0.47" ls 10,\
"../OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/0.46) w lp title "{/Times-Italic T} = 0.46" ls 11,\
"../OUTPUT/T0.45/N1080/noisecorr_NVT_combine_M5.txt" using 1:($2/0.45) w lp title "{/Times-Italic T} = 0.45" ls 12,\
0 ls 0 notitle


# Only diagonal correlation function, without normalizing, and including factor 1/T
reset
set output "./FIGURES/CdonT_unnorm.eps"
set xlabel "time [LJ units]"
set ylabel "{/Times-Italic C}_d(t)/T" offset 2.5,0
set key invert font ",27" samplen 2 maxrows 5
set tics font ",24"
set xtics format "10^{%T}"
load '~/.gnuplotting/gnuplot-palettes-master/rdbu.pal'
set logs x
plot   "../OUTPUT/T5.0/N1080/Cd_NVT.txt" using 1:($2/5.0) w lp title "{/Times-Italic T} = 5.0  " ls 1,\
"../OUTPUT/T2.0/N1080/Cd_NVT.txt" using 1:($2/2.0) w lp title "{/Times-Italic T} = 2.0  " ls 2,\
"../OUTPUT/T1.0/N1080/Cd_NVT.txt" using 1:($2/1.0) w lp title "{/Times-Italic T} = 1.0  " ls 3,\
"../OUTPUT/T0.8/N1080/Cd_NVT.txt" using 1:($2/0.8) w lp title "{/Times-Italic T} = 0.8  " ls 4,\
"../OUTPUT/T0.7/N1080/Cd_NVT.txt" using 1:($2/0.7) w lp title "{/Times-Italic T} = 0.7  " ls 5,\
"../OUTPUT/T0.6/N1080/Cd_NVT.txt" using 1:($2/0.6) w lp title "{/Times-Italic T} = 0.6  " ls 6,\
"../OUTPUT/T0.55/N1080/Cd_NVT.txt" using 1:($2/0.55) w lp title "{/Times-Italic T} = 0.55" ls 7,\
"../OUTPUT/T0.52/N1080/Cd_NVT.txt" using 1:($2/0.52) w lp title "{/Times-Italic T} = 0.52" ls 8,\
"../OUTPUT/T0.49/N1080/Cd_NVT.txt" using 1:($2/0.49) w lp title "{/Times-Italic T} = 0.49" ls 9,\
"../OUTPUT/T0.47/N1080/Cd_NVT.txt" using 1:($2/0.47) w lp title "{/Times-Italic T} = 0.47" ls 10,\
"../OUTPUT/T0.46/N1080/Cd_NVT.txt" using 1:($2/0.46) w lp title "{/Times-Italic T} = 0.46" ls 11,\
"../OUTPUT/T0.45/N1080/Cd_NVT.txt" using 1:($2/0.45) w lp title "{/Times-Italic T} = 0.45" ls 12,\
0 ls 0 notitle
set out





#####################################################
#                                                   #
# Two panels, comparing K(t) with Cd(t) and kMCT(t) #
#                                                   #
#####################################################
reset

set term pdfcairo enh font "Times-New-Roman,20" size 6,2
set out "./FIGURES/K-Cd-Kmct".kmax.".pdf"
#set term post enh c eps font "Times-Roman, 24" size 6,2
#set out "./FIGURES/K-Cd-Kmct.eps"
set logs x
load 'scl.pal'
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
normCd045=`awk '($1==0){print $2}' ../OUTPUT/T0.45/N1080/Cd_NVT.txt`

set multiplot layout 1,3

set border 0
set ylabel "Memory ({/Times-Italic t})" offset 19,1
unset key
unset xlabel
unset tics
set yr[-.1:1.05]
plot 1/0

marg=2.3
set lmargin marg
set tmargin .2
set rmargin 0
set bmargin 2.5
set xlabel "{/Times-Italic t}" offset 0,.8
unset key
#set key invert font ",12" samplen 0.5 maxrows 6
unset ylabel
set tics font "Times-New-Roman,12" nomirror
set xtics format "10^{%T}" offset -0.2,0.4
set ytics offset 0.5,0
set border
#set label "{/Times-Italic K}({/Times-Italic t})" at 0.003,0.1
set label "{/Times-Italic C}_{/Times-New-Roman d}({/Times-Italic t})" at 2000,0.2
set yr[-.1:1.05]
plot"../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 5.0  " ls 21 pt 1 ps .5 dt 1,\
"../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 2.0  " ls 22 pt 1 ps .5 dt 1,\
"../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 1.0  " ls 23 pt 1 ps .5 dt 1,\
"../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 0.8  " ls 24 pt 1 ps .5 dt 1,\
"../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 0.7  " ls 25 pt 1 ps .5 dt 1,\
"../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w lp notitle "{/Times-Italic T} = 0.6  " ls 26 pt 1 ps .5 dt 1,\
"../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.55" ls 27 pt 1 ps .5 dt 1,\
"../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.52" ls 28 pt 1 ps .5 dt 1,\
"../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.49" ls 29 pt 1 ps .5 dt 1,\
"../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.47" ls 30 pt 1 ps .5 dt 1,\
"../OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.46" ls 31 pt 1 ps .5 dt 1,\
"../OUTPUT/T0.45/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w lp notitle "{/Times-Italic T} = 0.45" ls 32 pt 1 ps .5 dt 1,\
0 ls 0 notitle,\
"../OUTPUT/T5.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd50):($3/normCd50) w lp  title "{/Times-Italic T} = 5.0  " ls 1 pt 2 ps .3 dt 1,\
"../OUTPUT/T2.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd20):($3/normCd20) w lp title "{/Times-Italic T} = 2.0  " ls 2 pt 2 ps .3 dt 1,\
"../OUTPUT/T1.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd10):($3/normCd10) w lp title "{/Times-Italic T} = 1.0  " ls 3 pt 2 ps .3 dt 1,\
"../OUTPUT/T0.8/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd08):($3/normCd08) w lp title "{/Times-Italic T} = 0.8  " ls 4 pt 2 ps .3 dt 1,\
"../OUTPUT/T0.7/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd07):($3/normCd07) w lp title "{/Times-Italic T} = 0.7  " ls 5 pt 2 ps .3 dt 1,\
"../OUTPUT/T0.6/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd06):($3/normCd06) w lp title "{/Times-Italic T} = 0.6  " ls 6 pt 2 ps .3 dt 1,\
"../OUTPUT/T0.55/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd055):($3/normCd055) w lp title "{/Times-Italic T} = 0.55" ls 7 pt 2 ps .3 dt 1,\
"../OUTPUT/T0.52/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd052):($3/normCd052) w lp title "{/Times-Italic T} = 0.52" ls 8 pt 2 ps .3 dt 1,\
"../OUTPUT/T0.49/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd049):($3/normCd049) w lp title "{/Times-Italic T} = 0.49" ls 9 pt 2 ps .3 dt 1,\
"../OUTPUT/T0.47/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd047):($3/normCd047) w lp title "{/Times-Italic T} = 0.47" ls 10 pt 2 ps .3 dt 1,\
"../OUTPUT/T0.46/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd046):($3/normCd046) w lp title "{/Times-Italic T} = 0.46" ls 11 pt 2 ps .3 dt 1,\
"../OUTPUT/T0.45/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd045):($3/normCd045) w lp title "{/Times-Italic T} = 0.45" ls 12 pt 2 ps .3 dt 1

set lmargin 0
set tmargin .2
set rmargin marg
set bmargin 2.5
unset ytics
set y2tics format "10^{%T}" mirror offset -0.7,0.2
set xtics format "10^{%T}" offset 0,0.4
unset label
unset key
set key invert font ",12" samplen 0.5 maxrows 6 top right at 1e4, 1.06
#set key invert font ",12" samplen 0.8 maxrows 6 at 1e-3,1.0
#set label "{/Times-Italic K}({/Times-Italic t})" at 0.003,0.1
set logs y2
set label "{/Times-Italic K}_{/Times-New-Roman MCT} ({/Times-Italic t})" at 10,0.4
set y2r[1e4:1.5e7]
plot[:1e4]"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 5.0  " ls 1  pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T2.0/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 2.0  " ls 2  pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T1.0/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 1.0  " ls 3  pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.8/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 0.8  " ls 4  pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.7/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 0.7  " ls 5  pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.6/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 0.6  " ls 6  pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.55/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 0.55" ls 7  pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.52/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 0.52" ls 8  pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.49/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 0.49" ls 9  pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.47/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 0.47" ls 10 pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 0.46" ls 11 pt 10 ps .5 dt 1,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.45/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp axis x1y2 t "{/Times-Italic T} = 0.45" ls 12 pt 10 ps .5 dt 1

width=0.22
height=0.37
set lmargin at screen 0.4460
set rmargin at screen 0.4460+width
set bmargin at screen 0.61
set tmargin at screen 0.61+height
unset xlabel
unset ylabel
unset logs y
set tics font ",12"
unset y2tics
set ytics format "%g" 0,.5,1 offset 0.7,-0.1
set xtics offset 0,0.6
unset label
set label "{/Times-Italic K}({/Times-Italic t})" at 5,0.7
set yr[-.1:1.05]
plot[:2e2]"../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w l notitle "{/Times-Italic T} = 5.0  " ls 1  dt 1,\
"../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w l notitle "{/Times-Italic T} = 2.0  " ls 2  dt 1,\
"../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w l notitle "{/Times-Italic T} = 1.0  " ls 3  dt 1,\
"../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w l notitle "{/Times-Italic T} = 0.8  " ls 4  dt 1,\
"../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w l notitle "{/Times-Italic T} = 0.7  " ls 5  dt 1,\
"../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M5.txt"  using 1:3 w l notitle "{/Times-Italic T} = 0.6  " ls 6  dt 1,\
"../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w l notitle "{/Times-Italic T} = 0.55" ls 7  dt 1,\
"../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w l notitle "{/Times-Italic T} = 0.52" ls 8  dt 1,\
"../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w l notitle "{/Times-Italic T} = 0.49" ls 9  dt 1,\
"../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w l notitle "{/Times-Italic T} = 0.47" ls 10  dt 1,\
"../OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w l notitle "{/Times-Italic T} = 0.46" ls 11  dt 1,\
"../OUTPUT/T0.45/N1080/noisecorr_NVT_combine_M5.txt" using 1:3 w l notitle "{/Times-Italic T} = 0.45" ls 12  dt 1,\
0 ls 0 notitle


unset multiplot


quit

