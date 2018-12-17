#!/usr/bin/gnuplot

#
# All temperatures together
#
set xlabel "t"
set ylabel "F_k(t)"
set xtics format "10^{%T}"
set logs x
set key bottom left invert
set title "N = 65"
plot [:][0:1] 0 linecolor -1 notitle, exp(-1) lc 0 lt 3 notitle,\
     "../OUTPUT/T10.0/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 10.0",\
     "../OUTPUT/T10.0/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 10.0",\
     "../OUTPUT/T5.0/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 5.0",\
     "../OUTPUT/T5.0/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 5.0",\
     "../OUTPUT/T2.0/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 2.0",\
     "../OUTPUT/T2.0/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 2.0",\
     "../OUTPUT/T1.0/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 1.0",\
     "../OUTPUT/T1.0/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 1.0",\
     "../OUTPUT/T0.8/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 0.8",\
     "../OUTPUT/T0.8/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 0.8",\
     "../OUTPUT/T0.7/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 0.7",\
     "../OUTPUT/T0.7/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 0.7",\
     "../OUTPUT/T0.6/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 0.6",\
     "../OUTPUT/T0.6/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 0.6",\
     "../OUTPUT/T0.55/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 0.55",\
     "../OUTPUT/T0.55/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 0.55",\
     "../OUTPUT/T0.52/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 0.52",\
     "../OUTPUT/T0.52/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 0.52",\
     "../OUTPUT/T0.5/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 0.5",\
     "../OUTPUT/T0.5/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 0.5",\
     "../OUTPUT/T0.49/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 0.49",\
     "../OUTPUT/T0.49/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 0.49",\
     "../OUTPUT/T0.48/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 0.48",\
     "../OUTPUT/T0.48/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 0.48",\
     "../OUTPUT/T0.47/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "T = 0.47",\
     "../OUTPUT/T0.47/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "T = 0.47"


set title "N = 1080"
plot [:][0:1] 0 linecolor -1 notitle, exp(-1) lc 0 lt 3 notitle,\
     "../OUTPUT/T5.0/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "T = 5.0",\
     "../OUTPUT/T5.0/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "T = 5.0",\
     "../OUTPUT/T2.0/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "T = 2.0",\
     "../OUTPUT/T2.0/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "T = 2.0",\
     "../OUTPUT/T1.0/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "T = 1.0",\
     "../OUTPUT/T1.0/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "T = 1.0",\
     "../OUTPUT/T0.8/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "T = 0.8",\
     "../OUTPUT/T0.8/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "T = 0.8",\
     "../OUTPUT/T0.7/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "T = 0.7",\
     "../OUTPUT/T0.7/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "T = 0.7",\
     "../OUTPUT/T0.6/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "T = 0.6",\
     "../OUTPUT/T0.6/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "T = 0.6",\
     "../OUTPUT/T0.55/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "T = 0.55",\
     "../OUTPUT/T0.55/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "T = 0.55"




#
# T=5.0
#
set title "T=5.0, N=1080"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T5.0/N1080/S".i."/Fkt_ifr0_xplor.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T5.0/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T5.0/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T5.0/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "after gap"
set title "T=5.0, N=65"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T5.0/N65/shift/S".i."/Fkt_ifr0_shift.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T5.0/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T5.0/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T5.0/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "after gap"


#
# T=2.0
#
set title "T=2.0, N=1080"
set xlabel "t"
set ylabel "F_k(t)"
p for [i=0:9] "../OUTPUT/T2.0/N1080/S".i."/Fkt_ifr0_xplor.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T2.0/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T2.0/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T2.0/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "after gap"
set title "T=2.0, N=65"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T2.0/N65/shift/S".i."/Fkt_ifr0_shift.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T2.0/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T2.0/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T2.0/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "after gap"

#
# T=1.0
#
set title "T=1.0, N=1080"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T1.0/N1080/S".i."/Fkt_ifr0_xplor.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T1.0/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T1.0/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T1.0/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "after gap"
set title "T=1.0, N=65"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T1.0/N65/shift/S".i."/Fkt_ifr0_shift.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T1.0/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T1.0/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T1.0/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "after gap"

#
# T=0.8
#
set title "T=0.8, N=1080"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T0.8/N1080/S".i."/Fkt_ifr0_xplor.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.8/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T0.8/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.8/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "after gap"
set title "T=0.8, N=65"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T0.8/N65/shift/S".i."/Fkt_ifr0_shift.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.8/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T0.8/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.8/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "after gap"

#
# T=0.7
#
set title "T=0.7"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T0.7/N1080/S".i."/Fkt_ifr0_xplor.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.7/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T0.7/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.7/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "after gap"
set title "T=0.7, N=65"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T0.7/N65/shift/S".i."/Fkt_ifr0_shift.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.7/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T0.7/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.7/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "after gap"



#
# T=0.6
#
## N=1080 ##
set title "T=0.6"
set xlabel "t"
set ylabel "F_k(t)"
set key bottom left
set logs x
set xtics format "10^{%T}"
p for [i=0:9] "../OUTPUT/T0.6/N1080/S".i."/Fkt_ifr0_xplor.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.6/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T0.6/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.6/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "after gap"
## N=  65 ##
set title "T=0.6, N=65"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T0.6/N65/shift/S".i."/Fkt_ifr0_shift.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.6/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T0.6/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.6/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "after gap"
# Compare NVE and NVT xplor
set logs x
p [][0:1]"../OUTPUT/T0.6/N65/xplor/Fkt_ifr0_xplor_NVT.txt" u 1:2:3 w errorl t "NVT First measurement",\
"../OUTPUT/T0.6/N65/xplor/Fkt_aftergap_xplor_NVT.txt" u 1:2:3 w errorl t "NVT after gap",\
"../OUTPUT/T0.6/N65/xplor/Fkt_ifr0_xplor_NVE.txt" u 1:2:3 w errorl t "NVE First measurement",\
"../OUTPUT/T0.6/N65/xplor/Fkt_aftergap_xplor_NVE.txt" u 1:2:3 w errorl t "NVE after gap"
# Compare NVE and NVT shift
set logs x
p [][0:1]"../OUTPUT/T0.6/N65/shift/Fkt_ifr0_shift_NVT.txt" u 1:2:3 w errorl t "NVT First measurement",\
"../OUTPUT/T0.6/N65/shift/Fkt_aftergap_shift_NVT.txt" u 1:2:3 w errorl t "NVT after gap",\
"../OUTPUT/T0.6/N65/shift/Fkt_ifr0_shift_NVE.txt" u 1:2:3 w errorl t "NVE First measurement",\
"../OUTPUT/T0.6/N65/shift/Fkt_aftergap_shift_NVE.txt" u 1:2:3 w errorl t "NVE after gap"





#
# T=0.55
#
set title "T=0.55"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T0.55/N1080/S".i."/Fkt_ifr0_xplor.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.55/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T0.55/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.55/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "after gap"
set title "T=0.55, N=65"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T0.55/N65/shift/S".i."/Fkt_ifr0_shift.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.55/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T0.55/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.55/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "after gap"

#
# T=0.52
#
set title "T=0.52"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T0.52/N1080/S".i."/Fkt_ifr0_xplor.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.52/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T0.52/N1080/Fkt_ifr0_xplor.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.52/N1080/Fkt_aftergap_xplor.txt" u 1:2:3 w errorl t "after gap"
set title "T=0.52, N=65"
set xlabel "t"
set ylabel "F_k(t)"
set logs x
p for [i=0:9] "../OUTPUT/T0.52/N65/shift/S".i."/Fkt_ifr0_shift.txt" u 1:2 w lp t "S".i
rep "../OUTPUT/T0.52/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl lw 4 lc 0 t "Average"
rep 0 lc -1
p "../OUTPUT/T0.52/N65/shift/Fkt_ifr0_shift.txt" u 1:2:3 w errorl t "First measurement",\
"../OUTPUT/T0.52/N65/shift/Fkt_aftergap_shift.txt" u 1:2:3 w errorl t "after gap"






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
