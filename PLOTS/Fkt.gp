#!/usr/bin/gnuplot

#
# All temperatures together
#
set xlabel "t"
set ylabel "F_k(t)"
set xtics format "10^{%T}"
set logs x
set key bottom left invert
set title "N = 65, thermalization check"

plot [:2e4][0:1] 0 linecolor -1 notitle, exp(-1) lc 0 lt 3 notitle,\
     "../OUTPUT/T10.0/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 10.0",\
     "../OUTPUT/T10.0/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "T = 10.0",\
     "../OUTPUT/T5.0/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 5.0",\
     "../OUTPUT/T5.0/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "T = 5.0",\
     "../OUTPUT/T2.0/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 2.0",\
     "../OUTPUT/T2.0/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "T = 2.0",\
     "../OUTPUT/T1.0/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 1.0",\
     "../OUTPUT/T1.0/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "T = 1.0",\
     "../OUTPUT/T0.8/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.8",\
     "../OUTPUT/T0.8/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.8",\
     "../OUTPUT/T0.7/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.7",\
     "../OUTPUT/T0.7/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.7",\
     "../OUTPUT/T0.6/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.6",\
     "../OUTPUT/T0.6/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.6",\
     "../OUTPUT/T0.55/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.55",\
     "../OUTPUT/T0.55/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.55",\
     "../OUTPUT/T0.52/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.52",\
     "../OUTPUT/T0.52/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.52",\
     "../OUTPUT/T0.49/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.49",\
     "../OUTPUT/T0.49/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.49"

plot [:][0:1] 0 linecolor -1 notitle, exp(-1) lc 0 lt 3 notitle,\
     "../OUTPUT/T0.6/N65/shift/Fkt_ifr0_shift_NVE.txt" u 1:2:3 w errorl t "T = 0.6",\
     "../OUTPUT/T0.49/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.49",\
     "../OUTPUT/T0.46/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "T = 0.46"


#
# T=0.6
#
set xlabel "t"
set ylabel "F_k(t)"
set logs x
set grid
set title "T=0.6, N=65, compare samples before gap"
p for [i=0:49] "../OUTPUT/T0.6/N65/shift/S".i."/Fkt_ifr0_shift_NVE.txt" u 1:2 w lp t "S".i,\
"../OUTPUT/T0.6/N65/shift/Fkt_ifr0_shift_NVE.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
set title "T=0.6, N=65, compare samples after gap"
p for [i=0:49] "../OUTPUT/T0.6/N65/shift/S".i."/Fkt_aftergap_shift_NVT.txt" u 1:2 w lp t "S".i,\
"../OUTPUT/T0.6/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
set title "T=0.6, N=65, check aging"
p "../OUTPUT/T0.6/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.6/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "after gap"


set xlabel "t"
set ylabel "F_k(t)"
set logs x
set grid
set title "T=0.49, N=65, compare samples before gap"
p for [i=0:49] "../OUTPUT/T0.49/N65/shift/S".i."/Fkt_ifr0_shift_NVT.txt" u 1:2 w lp t "S".i,\
"../OUTPUT/T0.49/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
set title "T=0.49, N=65, compare samples after gap"
p for [i=0:49] "../OUTPUT/T0.49/N65/shift/S".i."/Fkt_aftergap_shift_NVT.txt" u 1:2 w lp t "S".i,\
"../OUTPUT/T0.49/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
set title "T=0.49, N=65, check aging"
p "../OUTPUT/T0.49/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.49/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "after gap"


