#!/usr/bin/env gnuplot
reset
set term post enh c eps font "Times-Roman,32"

########################################
# Show noise autocorrelation functions #
########################################

set output "./FIGURES/K-Cd.eps"
set lmargin 6
set tmargin .5
set rmargin .5
set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation" offset 2.5,0
set key invert font ",27" samplen 2 maxrows 6
set tics font ",24"
set xtics format "10^{%T}"
set label "{/Times-Italic K}({/Times-Italic t})" at 0.005,0.1
set label "{/Times-Italic C}_d({/Times-Italic t})" at 100,0.3
load '~/.gnuplotting/gnuplot-palettes-master/rdbu.pal'
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
plot[][-.1:1.1]"../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "{/Times-Italic T} = 5.0  " ls 1 pt 1,\
"../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "{/Times-Italic T} = 2.0  " ls 2 pt 1,\
"../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "{/Times-Italic T} = 1.0  " ls 3 pt 1,\
"../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "{/Times-Italic T} = 0.8  " ls 4 pt 1,\
"../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "{/Times-Italic T} = 0.7  " ls 5 pt 1,\
"../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M3.txt"  using 1:3 w lp title "{/Times-Italic T} = 0.6  " ls 6 pt 1,\
"../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "{/Times-Italic T} = 0.55" ls 7 pt 1,\
"../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "{/Times-Italic T} = 0.52" ls 8 pt 1,\
"../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "{/Times-Italic T} = 0.49" ls 9 pt 1,\
"../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "{/Times-Italic T} = 0.47" ls 10 pt 1,\
"../OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M3.txt" using 1:3 w lp title "{/Times-Italic T} = 0.46" ls 11 pt 1,\
0 ls 0 notitle,\
"../OUTPUT/T5.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd50):($3/normCd50) w errorl notitle ls 1 pt 2,\
"../OUTPUT/T2.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd20):($3/normCd20) w errorl notitle ls 2 pt 2,\
"../OUTPUT/T1.0/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd10):($3/normCd10) w errorl notitle ls 3 pt 2,\
"../OUTPUT/T0.8/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd08):($3/normCd08) w errorl notitle ls 4 pt 2,\
"../OUTPUT/T0.7/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd07):($3/normCd07) w errorl notitle ls 5 pt 2,\
"../OUTPUT/T0.6/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd06):($3/normCd06) w errorl notitle ls 6 pt 2,\
"../OUTPUT/T0.55/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd055):($3/normCd055) w errorl notitle ls 7 pt 2,\
"../OUTPUT/T0.52/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd052):($3/normCd052) w errorl notitle ls 8 pt 2,\
"../OUTPUT/T0.49/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd049):($3/normCd049) w errorl notitle ls 9 pt 2,\
"../OUTPUT/T0.47/N1080/Cd_NVT.txt" using ($1*0.0025):($2/normCd047):($3/normCd047) w errorl notitle ls 10 pt 2
set out
unset label


# Only noise correlation function, without normalizing, and including factor 1/T
reset
set output "./FIGURES/KonT_unnorm.eps"
set xlabel "time [LJ units]"
set ylabel "Autocorrelation" offset 2.5,0
set key invert font ",27" samplen 2 maxrows 5
set tics font ",24"
set xtics format "10^{%T}"
load '~/.gnuplotting/gnuplot-palettes-master/rdbu.pal'
set logs x
plot   "../OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M3.txt" using 1:($2/5.0) w lp title "{/Times-Italic T} = 5.0  " ls 1,\
"../OUTPUT/T2.0/N1080/noisecorr_NVT_combine_M3.txt" using 1:($2/2.0) w lp title "{/Times-Italic T} = 2.0  " ls 2,\
"../OUTPUT/T1.0/N1080/noisecorr_NVT_combine_M3.txt" using 1:($2/1.0) w lp title "{/Times-Italic T} = 1.0  " ls 3,\
"../OUTPUT/T0.8/N1080/noisecorr_NVT_combine_M3.txt" using 1:($2/0.8) w lp title "{/Times-Italic T} = 0.8  " ls 4,\
"../OUTPUT/T0.7/N1080/noisecorr_NVT_combine_M3.txt" using 1:($2/0.7) w lp title "{/Times-Italic T} = 0.7  " ls 5,\
"../OUTPUT/T0.6/N1080/noisecorr_NVT_combine_M3.txt" using 1:($2/0.6) w lp title "{/Times-Italic T} = 0.6  " ls 6,\
"../OUTPUT/T0.55/N1080/noisecorr_NVT_combine_M3.txt" using 1:($2/0.55) w lp title "{/Times-Italic T} = 0.55" ls 7,\
"../OUTPUT/T0.52/N1080/noisecorr_NVT_combine_M3.txt" using 1:($2/0.52) w lp title "{/Times-Italic T} = 0.52" ls 8,\
"../OUTPUT/T0.49/N1080/noisecorr_NVT_combine_M3.txt" using 1:($2/0.49) w lp title "{/Times-Italic T} = 0.49" ls 9,\
"../OUTPUT/T0.47/N1080/noisecorr_NVT_combine_M3.txt" using 1:($2/0.47) w lp title "{/Times-Italic T} = 0.47" ls 10,\
0 ls 0 notitle






