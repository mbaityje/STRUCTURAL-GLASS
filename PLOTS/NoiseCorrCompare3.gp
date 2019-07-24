#!/usr/bin/env gnuplot
#
# Compare the K(t) with the vertex method to the other two #
#                                                          #
############################################################

###########
# T = 5.0 #
###########
set title "T=5.0"
unset logs
set xr[:.1]
set xlabel "t"
set ylabel "K(t)"
p "~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/Kvertex_NVT.txt" u ($1):3 w lp t "MCT-old" lw 2,\
"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/KvertexFuchs_NVT.txt" u ($1):3 w lp t "MCT-new" lw 2,\
"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/KvertexFuchs_NVT_kmax25.5.txt" u ($1):3 w lp t "MCT-new kmax" lw 2,\
"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/noisecorr_NVT.txt" u 1:3 w lp t "Exact",\
"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/Cd_NVT.txt" u ($1*0.0025):($2/1930.4532968941):($3/1930.4532968941) w errorl t "Diag"
set logs x
set xr[*:*]
rep



#Compare at time t=0
set title "T=5.0"
set xlabel "t"
set ylabel "K(t)"
unset logs x
p[:.01][:] "~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/Kvertex_NVT.txt" u ($1):($2) w lp t "MCT-old",\
"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/KvertexFuchs_NVT.txt" u ($1):($2) w lp t "MCT-new",\
"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/noisecorr_NVT.txt" u 1:2 w lp t "Exact"




############
# T = 0.46 #
############
set title "T=0.46"
unset logs
set xr[:.1]
set xlabel "t"
set ylabel "K(t)"
p "~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/Kvertex_NVT.txt" u ($1):3 w lp t "MCT-old" lw 2,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/KvertexFuchs_NVT.txt" u ($1):3 w lp t "MCT-new" lw 2,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/KvertexFuchs_NVT_kmax25.5.txt" u ($1):3 w lp t "MCT-new kmax" lw 2,\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/noisecorr_NVT.txt" u 1:3 w lp t "Exact",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/Cd_NVT.txt" u ($1*0.0025):($2/667.41270962594):($3/667.41270962594) w errorl t "Diag"
set logs x
set xr[*:*]
rep

#Compare at time t=0
set title "T=0.46"
set xlabel "t"
set ylabel "K(t)"
unset logs x
p [:0.01]"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/Kvertex_NVT.txt" u ($1):($2) w lp t "MCT-old",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/KvertexFuchs_NVT.txt" u ($1):($2) w lp t "MCT-new",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/noisecorr_NVT.txt" u 1:($2/0.46) w lp t "Exact"
#rep "~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/CFF_NVT.txt" u 1:2 w lp t "CFF"


#########################
# T=5.0, 0.46 multiplot #
#########################
reset
set term post enh c eps size 3,1.1 font "Times-Roman,14" dashed
set out "./FIGURES/Kmct.eps"
set multiplot layout 1,2
set yr[-.05:1]
set logs x
set tics font ",10" nomirror
set xtics format "10^{%T}" offset 0,0.45
set xlabel "{/Times-Italic t}" offset 0,1 font ",16"
set tmargin 0.5
set bmargin 2
side=4.5

set ylabel "Memory ({/Times-Italic t})" offset 3,0
set ytics offset 0.5,0
set lmargin side
set rmargin 0
set key samplen 1 at 6,0.8 font ",14" spacing 1.2
set label "{/Times-Italic T} = 5.0" at .02,.85 font ",14"
set label "(a)" at 4e-3,.1 font ",14"
p "~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/Cd_NVT.txt" u ($1*0.0025):($2/1930.4532968941):($3/1930.4532968941) w errorl ps 0.5 lt rgb "dark-green" t "{/Times-Italic C}_d({/Times-Italic t})",\
"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ps 0.4 lt rgb "dark-violet" t "{/Times-Italic K}({/Times-Italic t})",\
"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/KvertexFuchs_NVT_kmax28.txt" u ($1):3 w lp lw 2 ps 0.5 lc rgb "dark-orange" lt 4 t "{/Times-Italic K}_{MCT}({/Times-Italic t})",\
0 lc -1 lt 0 not
#"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/KvertexFuchs_NVT_kmax25.5.txt" u ($1):3 w lp lw 2 ps 0.5 lc rgb "dark-khaki" lt 6 t "{/Times-Italic K}@^*_{MCT}({/Times-Italic t})",\
#"~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/Kvertex_NVT.txt" u ($1):3 w lp lw 2 ps 0.5 lc rgb "dark-salmon" lt 2 t "{/Times-Italic K}_{MCT}({/Times-Italic t})",\

set lmargin 0
set rmargin side
set ytics format ""
unset label
set label "{/Times-Italic T} = 0.46" at .1,.85 font ",14"
set label "(b)" at 8e-3,.1 font ",14"
unset ylabel
unset key
p "~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/Cd_NVT.txt" u ($1*0.0025):($2/667.41270962594):($3/667.41270962594) w errorl ps 0.5 lt rgb "dark-green" t "{/Times-Italic C}_d({/Times-Italic t})",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/noisecorr_NVT_combine_M5.txt" u 1:3 w lp ps 0.4 lt rgb 'dark-violet' t "{/Times-Italic K}({/Times-Italic t})",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/KvertexFuchs_NVT_kmax28.txt" u ($1):3 w lp lw 2 ps 0.5 lc rgb "dark-orange" lt 4 t "{/Times-Italic K}_{MCT}({/Times-Italic t})",\
0 lc -1 lt 0 not
#"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/KvertexFuchs_NVT_kmax25.5.txt" u ($1):3 w lp lw 2 ps 0.5 lc rgb "dark-khaki" lt 6 t "{/Times-Italic K}@^*_{MCT}({/Times-Italic t})",\
#"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/Kvertex_NVT.txt" u ($1):3 w lp lw 2 ps 0.5 lc rgb "dark-salmon" lt 2 t "{/Times-Italic K}_{MCT}({/Times-Italic t})",\
set lmargin at screen 0.7
set rmargin at screen 0.905
set bmargin at screen 0.51
set tmargin at screen 0.92
unset label
set xtics format ""
set ytics ("10^{-4}" 1e-4, "10^{-2}" 1e-2, "1" 1)
set logs
set yr[1e-2:1]
unset xlabel
set lmargin at screen 0.7
#replot
unset multiplot


