set key top left
set title "N=65"
set xlabel "Temperature"
set ylabel "Energy"
p "SMALL-DATA/E_T.txt" u 2:3 w lp ps 1 pt 6 t "E(T)",\
-380 t "E_{th}" ls 0,\
115*x-415 t "linear extrapolation"




pause -1