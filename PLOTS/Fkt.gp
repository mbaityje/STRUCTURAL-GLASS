#!/usr/bin/gnuplot

#
# All temperatures together
#
set xlabel "t"
set ylabel "F_k(t)"
set xtics format "10^{%T}"
set logs x
set key bottom left invert
plot [:1e8][0:1] 0 linecolor -1 notitle, exp(-1) lc 0 lt 3 notitle,\
     "../OUTPUT/T10.0/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "T = 10.0",\
     "../OUTPUT/T10.0/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "T = 10.0",\
     "../OUTPUT/T2.0/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "T = 2.0",\
     "../OUTPUT/T2.0/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "T = 2.0",\
     "../OUTPUT/T0.6/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "T = 0.6",\
     "../OUTPUT/T0.6/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "T = 0.6",\
     "../OUTPUT/T0.49/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "T = 0.49",\
     "../OUTPUT/T0.49/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "T = 0.49",\
     "../OUTPUT/T0.466/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "T = 0.466",\
     "../OUTPUT/T0.466/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "T = 0.466",\
     "../OUTPUT/T0.44/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "T = 0.44",\
     "../OUTPUT/T0.44/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "T = 0.44",\
     "../OUTPUT/T0.43/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "T = 0.43",\
     "../OUTPUT/T0.43/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "T = 0.43",\
     "../OUTPUT/T0.42/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "T = 0.42",\
     "../OUTPUT/T0.42/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "T = 0.42"




#
# T=10.0
#
set title "T=10.0"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T10.0/N65/S".i."/Fkt_ifr0.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T10.0/N65/Fkt_ifr0.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1

p "../OUTPUT/T10.0/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T10.0/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "after gap"


#
# T=2.0
#
set title "T=2.0"
set xlabel "t"
set ylabel "F_k(t)"
p for [i=0:9] "../OUTPUT/T2.0/N65/S".i."/Fkt_ifr0.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T2.0/N65/Fkt_ifr0.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1

p "../OUTPUT/T2.0/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T2.0/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "after gap"




#
# T=0.6
#
set title "T=0.6"
set xlabel "t"
set ylabel "F_k(t)"
set key bottom left
set logs x
set xtics format "10^{%T}"
#Before gap
p [][0:1] for [i=0:39]"../OUTPUT/T0.6/N65/S".i."/Fkt_ifr0_xplor.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.6/N65/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
#After gap
p [][0:1] for [i=0:9]"../OUTPUT/T0.6/N65/S".i."/Fkt_aftergap_xplor.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.6/N65/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
#Compare averages
p [][:1]"../OUTPUT/T0.6/N65/Fkt_ifr0_10samples_xplor.txt" u 1:2:3 w errorl t "First measurement (10 samples)",\
"../OUTPUT/T0.6/N65/Fkt_aftergap_10samples_xplor.txt" u 1:2:3 w errorl t "after gap (10 samples)"
rep "../OUTPUT/T0.6/N65/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "First measurement (all samples)",\
"../OUTPUT/T0.6/N65/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "after gap (all samples)"


#
# T=0.49
#
set title "T=0.49"
set xlabel "t"
set ylabel "F_k(t)"
set key bottom left
set logs x
set xtics format "10^{%T}"
#Before gap
p [][0:1] for [i=0:9]"../OUTPUT/T0.49/N65/S".i."/Fkt_ifr0.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.49/N65/Fkt_ifr0.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
#After gap
p [][0:1] for [i=0:9]"../OUTPUT/T0.49/N65/S".i."/Fkt_aftergap.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.49/N65/Fkt_aftergap.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
#Compare averages
p [][0:1]"../OUTPUT/T0.49/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.49/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "after gap"



#
# T=0.43
#
set title "T=0.43"
set xlabel "t"
set ylabel "F_k(t)"
p for [i=0:9] "../OUTPUT/T0.43/N65/S".i."/Fkt_ifr0.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.43/N65/Fkt_ifr0.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"

p [][0:1]"../OUTPUT/T0.43/N65/Fkt_ifr0.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.43/N65/Fkt_aftergap.txt" u 1:2:3 w errorl t "after gap"
