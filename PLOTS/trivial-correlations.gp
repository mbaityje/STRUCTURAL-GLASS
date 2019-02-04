#!/usr/bin/env gnuplot
set term post enh c eps font "Times-Roman,32"


#Self-intermediate scattering function
set out "./FIGURES/Fkt.eps"
set lmargin 6
set rmargin .5
set tmargin .5
set ylabel "{/Times-Italic F_k}({/Times-Italic t})" offset 2,0
set xlabel "{/Times-Italic t}" offset 0,0
set logs x
set xtics format "10^{%T}" font ",24"
set ytics 0,0.2 font ",24"
set key font ",24" samplen 1 invert at 0.0007,0.38
load '~/.gnuplotting/gnuplot-palettes-master/rdbu.pal'
p "../OUTPUT/T5.0/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 5.0  " ls 1,\
"../OUTPUT/T2.0/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 2.0  " ls 2,\
"../OUTPUT/T1.0/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 1.0  " ls 3,\
"../OUTPUT/T0.8/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 0.8  " ls 4,\
"../OUTPUT/T0.7/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 0.7  " ls 5,\
"../OUTPUT/T0.6/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 0.6  " ls 6,\
"../OUTPUT/T0.55/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 0.55" ls 7,\
"../OUTPUT/T0.55/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 0.55" ls 8,\
"../OUTPUT/T0.52/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 0.52" ls 9,\
"../OUTPUT/T0.49/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 0.49" ls 10,\
"../OUTPUT/T0.47/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 0.47" ls 11,\
"../OUTPUT/T0.46/N1080/Fkt_NVT.txt"  using 1:2:3 w errorl title "{/Times-Italic T} = 0.46" ls 12
reset


#Force-Force correlator

set out "./FIGURES/CFF.eps"
set lmargin 6
set rmargin .5
set tmargin .5
set ylabel "{/Times-Italic C}^{FF}({/Times-Italic t})/{/Times-Italic C}^{FF}(0)" offset 2,0
set xlabel "{/Times-Italic t}" offset 0,0
set logs x
set xtics format "10^{%T}" font ",24"
set ytics 0,0.2 font ",24"
set key top right font ",24" samplen 1 invert maxrows 6
load '~/.gnuplotting/gnuplot-palettes-master/rdbu.pal'
normCFF50=`awk '($1==0){print $2}' ../OUTPUT/T5.0/N1080/CFF_NVT.txt`
normCFF20=`awk '($1==0){print $2}' ../OUTPUT/T2.0/N1080/CFF_NVT.txt`
normCFF10=`awk '($1==0){print $2}' ../OUTPUT/T1.0/N1080/CFF_NVT.txt`
normCFF08=`awk '($1==0){print $2}' ../OUTPUT/T0.8/N1080/CFF_NVT.txt`
normCFF07=`awk '($1==0){print $2}' ../OUTPUT/T0.7/N1080/CFF_NVT.txt`
normCFF06=`awk '($1==0){print $2}' ../OUTPUT/T0.6/N1080/CFF_NVT.txt`
normCFF055=`awk '($1==0){print $2}' ../OUTPUT/T0.55/N1080/CFF_NVT.txt`
normCFF052=`awk '($1==0){print $2}' ../OUTPUT/T0.52/N1080/CFF_NVT.txt`
normCFF049=`awk '($1==0){print $2}' ../OUTPUT/T0.49/N1080/CFF_NVT.txt`
normCFF047=`awk '($1==0){print $2}' ../OUTPUT/T0.47/N1080/CFF_NVT.txt`
normCFF046=`awk '($1==0){print $2}' ../OUTPUT/T0.46/N1080/CFF_NVT.txt`
p "../OUTPUT/T5.0/N1080/CFF_NVT.txt"  using 1:($2/normCFF50):($3/normCFF50) w errorl title "{/Times-Italic T} = 5.0  " ls 1,\
"../OUTPUT/T2.0/N1080/CFF_NVT.txt"  using 1:($2/normCFF20):($3/normCFF20) w errorl title "{/Times-Italic T} = 2.0  " ls 2,\
"../OUTPUT/T1.0/N1080/CFF_NVT.txt"  using 1:($2/normCFF10):($3/normCFF10) w errorl title "{/Times-Italic T} = 1.0  " ls 3,\
"../OUTPUT/T0.8/N1080/CFF_NVT.txt"  using 1:($2/normCFF08):($3/normCFF08) w errorl title "{/Times-Italic T} = 0.8  " ls 4,\
"../OUTPUT/T0.7/N1080/CFF_NVT.txt"  using 1:($2/normCFF07):($3/normCFF07) w errorl title "{/Times-Italic T} = 0.7  " ls 5,\
"../OUTPUT/T0.6/N1080/CFF_NVT.txt"  using 1:($2/normCFF06):($3/normCFF06) w errorl title "{/Times-Italic T} = 0.6  " ls 6,\
"../OUTPUT/T0.55/N1080/CFF_NVT.txt"  using 1:($2/normCFF055):($3/normCFF055) w errorl title "{/Times-Italic T} = 0.55" ls 7,\
"../OUTPUT/T0.52/N1080/CFF_NVT.txt"  using 1:($2/normCFF052):($3/normCFF052) w errorl title "{/Times-Italic T} = 0.52" ls 8,\
"../OUTPUT/T0.49/N1080/CFF_NVT.txt"  using 1:($2/normCFF049):($3/normCFF049) w errorl title "{/Times-Italic T} = 0.49" ls 9,\
"../OUTPUT/T0.47/N1080/CFF_NVT.txt"  using 1:($2/normCFF047):($3/normCFF047) w errorl title "{/Times-Italic T} = 0.47" ls 10,\
"../OUTPUT/T0.46/N1080/CFF_NVT.txt"  using 1:($2/normCFF046):($3/normCFF046) w errorl title "{/Times-Italic T} = 0.46" ls 11



# Momentum-momentum correlator

set out "./FIGURES/CPP.eps"
set lmargin 6
set rmargin .5
set tmargin .5
set ylabel "{/Times-Italic C}^{PP}({/Times-Italic t})/{/Times-Italic C}^{PP}(0)" offset 2,0
set xlabel "{/Times-Italic t}" offset 0,0
set logs x
set xtics format "10^{%T}" font ",24"
set ytics 0,0.2 font ",24"
set key top right font ",24" samplen 1 invert maxrows 6
load '~/.gnuplotting/gnuplot-palettes-master/rdbu.pal'
normCPP50=`awk '($1==0){print $2}' ../OUTPUT/T5.0/N1080/CPP_NVT.txt`
normCPP20=`awk '($1==0){print $2}' ../OUTPUT/T2.0/N1080/CPP_NVT.txt`
normCPP10=`awk '($1==0){print $2}' ../OUTPUT/T1.0/N1080/CPP_NVT.txt`
normCPP08=`awk '($1==0){print $2}' ../OUTPUT/T0.8/N1080/CPP_NVT.txt`
normCPP07=`awk '($1==0){print $2}' ../OUTPUT/T0.7/N1080/CPP_NVT.txt`
normCPP06=`awk '($1==0){print $2}' ../OUTPUT/T0.6/N1080/CPP_NVT.txt`
normCPP055=`awk '($1==0){print $2}' ../OUTPUT/T0.55/N1080/CPP_NVT.txt`
normCPP052=`awk '($1==0){print $2}' ../OUTPUT/T0.52/N1080/CPP_NVT.txt`
normCPP049=`awk '($1==0){print $2}' ../OUTPUT/T0.49/N1080/CPP_NVT.txt`
normCPP047=`awk '($1==0){print $2}' ../OUTPUT/T0.47/N1080/CPP_NVT.txt`
normCPP046=`awk '($1==0){print $2}' ../OUTPUT/T0.46/N1080/CPP_NVT.txt`
p "../OUTPUT/T5.0/N1080/CPP_NVT.txt"  using 1:($2/normCPP50):($3/normCPP50) w errorl title "{/Times-Italic T} = 5.0  " ls 1,\
"../OUTPUT/T2.0/N1080/CPP_NVT.txt"  using 1:($2/normCPP20):($3/normCPP20) w errorl title "{/Times-Italic T} = 2.0  " ls 2,\
"../OUTPUT/T1.0/N1080/CPP_NVT.txt"  using 1:($2/normCPP10):($3/normCPP10) w errorl title "{/Times-Italic T} = 1.0  " ls 3,\
"../OUTPUT/T0.8/N1080/CPP_NVT.txt"  using 1:($2/normCPP08):($3/normCPP08) w errorl title "{/Times-Italic T} = 0.8  " ls 4,\
"../OUTPUT/T0.7/N1080/CPP_NVT.txt"  using 1:($2/normCPP07):($3/normCPP07) w errorl title "{/Times-Italic T} = 0.7  " ls 5,\
"../OUTPUT/T0.6/N1080/CPP_NVT.txt"  using 1:($2/normCPP06):($3/normCPP06) w errorl title "{/Times-Italic T} = 0.6  " ls 6,\
"../OUTPUT/T0.55/N1080/CPP_NVT.txt"  using 1:($2/normCPP055):($3/normCPP055) w errorl title "{/Times-Italic T} = 0.55" ls 7,\
"../OUTPUT/T0.52/N1080/CPP_NVT.txt"  using 1:($2/normCPP052):($3/normCPP052) w errorl title "{/Times-Italic T} = 0.52" ls 8,\
"../OUTPUT/T0.49/N1080/CPP_NVT.txt"  using 1:($2/normCPP049):($3/normCPP049) w errorl title "{/Times-Italic T} = 0.49" ls 9,\
"../OUTPUT/T0.47/N1080/CPP_NVT.txt"  using 1:($2/normCPP047):($3/normCPP047) w errorl title "{/Times-Italic T} = 0.47" ls 10,\
"../OUTPUT/T0.46/N1080/CPP_NVT.txt"  using 1:($2/normCPP046):($3/normCPP046) w errorl title "{/Times-Italic T} = 0.46" ls 11



# Force-momentum correlator

set out "./FIGURES/CFP.eps"
set lmargin 6
set rmargin .5
set tmargin .5
set ylabel "{/Times-Italic C}^{FP}({/Times-Italic t})" offset 2,0
set xlabel "{/Times-Italic t}" offset 0,0
set logs x
set xtics format "10^{%T}" font ",24"
set ytics auto font ",24"
set key top right font ",24" samplen 1 invert maxrows 6
load '~/.gnuplotting/gnuplot-palettes-master/rdbu.pal'
p "../OUTPUT/T5.0/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 5.0  " ls 1,\
"../OUTPUT/T2.0/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 2.0  " ls 2,\
"../OUTPUT/T1.0/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 1.0  " ls 3,\
"../OUTPUT/T0.8/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 0.8  " ls 4,\
"../OUTPUT/T0.7/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 0.7  " ls 5,\
"../OUTPUT/T0.6/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 0.6  " ls 6,\
"../OUTPUT/T0.55/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 0.55" ls 7,\
"../OUTPUT/T0.52/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 0.52" ls 8,\
"../OUTPUT/T0.49/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 0.49" ls 9,\
"../OUTPUT/T0.47/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 0.47" ls 10,\
"../OUTPUT/T0.46/N1080/CFP_NVT.txt"  using 1:($2):($3) w errorl title "{/Times-Italic T} = 0.46" ls 11

