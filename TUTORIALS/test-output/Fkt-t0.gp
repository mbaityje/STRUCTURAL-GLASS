reset
set xlabel "time step"
set ylabel "T"
set tics format "%g"
p "observables-test_NVT5000.txt" u 1:2, 1 t "T_{target}"

reset
set xlabel "time step"
set ylabel "energy"
set tics format "%g"
p "observables-test_NVT5000.txt" u 1:($4+$5)

reset
set xlabel "time step"
set ylabel "mean square displacement"
set key top left
set tics format "%g"
set logs
p "msd.txt" u ($1*5000):2 w l notitle,\
x>300000 ? 3e1*x**0.5 : 1/0 lw 4 t "x^{0.5}",\
x<300000 && x>80000 ? 1e-2*x : 1/0 lw 4 t "x"

reset
load '~/.gnuplotting/gnuplot-palettes-master/blues.pal'
set xlabel "time step (dt=0.001)"
set ylabel "F_k(t_0,t)"
set key bottom left
set tics format "%g"
set logs
p for [i=0:9]"< egrep ^".i." Fkt-t0.txt" u (($4-$3)*5000):5 ls i w l

