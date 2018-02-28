#!/usr/bin/gnuplot

set title "Force correlation, T=0.6"
set logs x
set xlabel "t"
set ylabel "C(t)/C(0)"
p "../OUTPUT/T0.6/N1080/S0/corrF.txt" u 2:5 w lp t "Full" lc 1
rep "../OUTPUT/T0.6/N1080/S1/corrF.txt" u 2:5 w lp t "Full" lc 1
rep "../OUTPUT/T0.6/N1080/S2/corrF.txt" u 2:5 w lp t "Full" lc 1
rep "../OUTPUT/T0.6/N1080/S3/corrF.txt" u 2:5 w lp t "Full" lc 1
rep "../OUTPUT/T0.6/N1080/S0/corrF.txt" u 2:6 w lp t "Diagonal" lc 2
rep "../OUTPUT/T0.6/N1080/S1/corrF.txt" u 2:6 w lp t "Diagonal" lc 2
rep "../OUTPUT/T0.6/N1080/S2/corrF.txt" u 2:6 w lp t "Diagonal" lc 2
rep "../OUTPUT/T0.6/N1080/S3/corrF.txt" u 2:6 w lp t "Diagonal" lc 2


set title "Force correlation, T=0.466"
set logs x
set xlabel "t"
set ylabel "C(t)/C(0)"
p "../OUTPUT/T0.466/N1080/S0/corrF.txt" u 2:5 w lp t "Full" lc 1
rep "../OUTPUT/T0.466/N1080/S1/corrF.txt" u 2:5 w lp t "Full" lc 1
rep "../OUTPUT/T0.466/N1080/S2/corrF.txt" u 2:5 w lp t "Full" lc 1
rep "../OUTPUT/T0.466/N1080/S3/corrF.txt" u 2:5 w lp t "Full" lc 1
rep "../OUTPUT/T0.466/N1080/S0/corrF.txt" u 2:6 w lp t "Diagonal" lc 2
rep "../OUTPUT/T0.466/N1080/S1/corrF.txt" u 2:6 w lp t "Diagonal" lc 2
rep "../OUTPUT/T0.466/N1080/S2/corrF.txt" u 2:6 w lp t "Diagonal" lc 2
rep "../OUTPUT/T0.466/N1080/S3/corrF.txt" u 2:6 w lp t "Diagonal" lc 2




