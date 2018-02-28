#!/usr/bin/gnuplot

#T=10
set logs x
set xlabel "t"
set title "T=10"
unset key

set ylabel "E"
p  "../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep0.txt" u 2:3 w lp
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep1.txt" u 2:3 w lp
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep2.txt" u 2:3 w lp
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep3.txt" u 2:3 w lp
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep4.txt" u 2:3 w lp
set ylabel "E_{IS}"
rep  "../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep0.txt" u 2:4 w lp
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep1.txt" u 2:4 w lp
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep2.txt" u 2:4 w lp
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep3.txt" u 2:4 w lp
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep4.txt" u 2:4 w lp
set ylabel "q, q_{IS}"
set key bottom left
p  "../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep0.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep1.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep2.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep3.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep4.txt" u 2:5 w lp lc 1 title "q"
rep  "../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep0.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep1.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep2.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep3.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T10.0/N65/S0/OverlapTrajectoryISirep4.txt" u 2:6 w lp lc 2 title "q_{IS}"




#T=2.0
set logs x
set xlabel "t"
set title "T=2.0"
unset key

set ylabel "E"
p  "../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep0.txt" u 2:3 w lp
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep1.txt" u 2:3 w lp
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep2.txt" u 2:3 w lp
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep3.txt" u 2:3 w lp
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep4.txt" u 2:3 w lp
set ylabel "E_{IS}"
p  "../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep0.txt" u 2:4 w lp
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep1.txt" u 2:4 w lp
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep2.txt" u 2:4 w lp
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep3.txt" u 2:4 w lp
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep4.txt" u 2:4 w lp
set ylabel "q, q_{IS}"
set key bottom left
p  "../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep0.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep1.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep2.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep3.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep4.txt" u 2:5 w lp lc 1 title "q"
rep  "../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep0.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep1.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep2.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep3.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T2.0/N65/S0/OverlapTrajectoryISirep4.txt" u 2:6 w lp lc 2 title "q_{IS}"






#T=0.43
set logs x
set xlabel "t"
set title "T=0.43"
unset key

set ylabel "E"
p  "../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep0.txt" u 2:3 w lp
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep1.txt" u 2:3 w lp
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep2.txt" u 2:3 w lp
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep3.txt" u 2:3 w lp
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep4.txt" u 2:3 w lp
set ylabel "E_{IS}"
p  "../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep0.txt" u 2:4 w lp
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep1.txt" u 2:4 w lp
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep2.txt" u 2:4 w lp
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep3.txt" u 2:4 w lp
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep4.txt" u 2:4 w lp
set ylabel "q, q_{IS}"
set key bottom left
p  "../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep0.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep1.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep2.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep3.txt" u 2:5 w lp lc 1 title "q"
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep4.txt" u 2:5 w lp lc 1 title "q"
rep  "../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep0.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep1.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep2.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep3.txt" u 2:6 w lp lc 2 title "q_{IS}"
rep"../OUTPUT/T0.43/N65/S0/OverlapTrajectoryISirep4.txt" u 2:6 w lp lc 2 title "q_{IS}"
