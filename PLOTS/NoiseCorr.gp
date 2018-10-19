#!/usr/bin/env gnuplot

set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation"
set key invert
load '~/.gnuplotting/gnuplot-palettes-master/moreland.pal'
load '~/.gnuplotting/gnuplot-palettes-master/rdbu.pal'
plot   "../OUTPUT/T5.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "T = 5.0" ls 1
replot "../OUTPUT/T2.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "T = 2.0" ls 2
replot "../OUTPUT/T1.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "T = 1.0" ls 3
replot "../OUTPUT/T0.6/N1080/noisecorr_NVE.txt" using 1:3 w lp title "T = 0.6" ls 4
replot 0 ls 0 notitle

#Without normalizing, and putting factor 1/T
plot   "../OUTPUT/T5.0/N1080/noisecorr_NVE.txt" using 1:($2/5.) w lp title "T = 5.0" ls 1
replot "../OUTPUT/T2.0/N1080/noisecorr_NVE.txt" using 1:($2/2.) w lp title "T = 2.0" ls 2
replot "../OUTPUT/T1.0/N1080/noisecorr_NVE.txt" using 1:($2/1.) w lp title "T = 1.0" ls 3
replot "../OUTPUT/T0.6/N1080/noisecorr_NVE.txt" using 1:($2/.6) w lp title "T = 0.6" ls 4
replot 0 ls 0 notitle


set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation"
load '~/.gnuplotting/gnuplot-palettes-master/paired.pal'
plot   "../OUTPUT/T5.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "K^{NVE}" ls 1
replot   "../OUTPUT/T5.0/N1080/noisecorrlinear_NVE.txt" using 1:3 w lp title "K_{l}^{NVE}" ls 3
replot "../OUTPUT/T5.0/N1080/noisecorrlinear-selfcon_NVE.txt" using 1:3 w lp title "K_{lsc}^{NVE}" ls 5
replot 0 ls 0 notitle


set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation"
load '~/.gnuplotting/gnuplot-palettes-master/paired.pal'
plot   "../OUTPUT/T2.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "K^{NVE}" ls 1
replot   "../OUTPUT/T2.0/N1080/noisecorrlinear_NVE.txt" using 1:3 w lp title "K_{l}^{NVE}" ls 3
replot "../OUTPUT/T2.0/N1080/noisecorrlinear-selfcon_NVE.txt" using 1:3 w lp title "K_{lsc}^{NVE}" ls 5
replot 0 ls 0 notitle


set logs x
set xlabel "time [LJ units]"
set ylabel "Autocorrelation"
load '~/.gnuplotting/gnuplot-palettes-master/paired.pal'
plot   "../OUTPUT/T1.0/N1080/noisecorr_NVE.txt" using 1:3 w lp title "K^{NVE}" ls 1
replot   "../OUTPUT/T1.0/N1080/noisecorrlinear_NVE.txt" using 1:3 w lp title "K_{l}^{NVE}" ls 3
replot "../OUTPUT/T1.0/N1080/noisecorrlinear-selfcon_NVE.txt" using 1:3 w lp title "K_{lsc}^{NVE}" ls 5
replot 0 ls 0 notitle


