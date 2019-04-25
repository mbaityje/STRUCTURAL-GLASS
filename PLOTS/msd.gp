#!/usr/bin/gnuplot

set term post enh c eps
set out "./FIGURES/msd.eps"
set logs
set key bottom right
set xlabel "{/Times-Italic t}"
set ylabel "Mean Square Displacement"
set tics format "10^{%T}"
p "../OUTPUT/T5.0/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 5.0",\
"../OUTPUT/T2.0/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 2.0",\
"../OUTPUT/T1.0/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 1.0",\
"../OUTPUT/T0.8/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 0.8",\
"../OUTPUT/T0.7/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 0.7",\
"../OUTPUT/T0.6/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 0.6",\
"../OUTPUT/T0.55/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 0.55",\
"../OUTPUT/T0.52/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 0.52",\
"../OUTPUT/T0.49/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 0.49",\
"../OUTPUT/T0.47/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 0.47",\
"../OUTPUT/T0.46/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 0.46",\
"../OUTPUT/T0.45/N1080/msd_NVT.txt" u 1:2:3 w errorl t "T = 0.45"
set out

#Find the diffusion coefficients through msd=D/6t at large t
reset
set out "./FIGURES/msdont.eps"
set logs
set key bottom right
set xlabel "{/Times-Italic t}"
set ylabel "{/Symbol D}^2 / {/Times-Italic t}"
plot "../OUTPUT/T5.0/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 5.0",\
"../OUTPUT/T2.0/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 2.0",\
"../OUTPUT/T1.0/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 1.0",\
"../OUTPUT/T0.8/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.8",\
"../OUTPUT/T0.7/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.7",\
"../OUTPUT/T0.6/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.6",\
"../OUTPUT/T0.55/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.55",\
"../OUTPUT/T0.52/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.52",\
"../OUTPUT/T0.49/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.49",\
"../OUTPUT/T0.47/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.47",\
"../OUTPUT/T0.46/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.46",\
"../OUTPUT/T0.45/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.45"

f50(x)=D50+a50/x
f20(x)=D20+a20/x
f10(x)=D10+a10/x
f08(x)=D08+a08/x
f07(x)=D07+a07/x
f06(x)=D06+a06/x
f055(x)=D055+a055/x
f052(x)=D052+a052/x
f049(x)=D049+a049/x
f047(x)=D047+a047/x
f046(x)=D046+a046/x
f045(x)=D045+a045/x
D50=0.6;     a50=1
D20=0.3;     a20=1
D10=0.07;    a10=1
D08=0.035;   a08=1
D07=0.016;   a07=1
D06=0.008;   a06=1
D055=0.004;  a55=1
D052=0.002;  a52=1
D049=0.001;  a49=1
D047=0.0005; a47=1
D046=0.0002; a46=1
D045=0.0001; a45=1

set fit errorvariables
unset logs
set logs x

fit [3:] f50(x) "../OUTPUT/T5.0/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D50,a50
plot [.1:][:]"../OUTPUT/T5.0/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 5.0",f50(x)

fit [9:] f20(x) "../OUTPUT/T2.0/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D20,a20
plot [1.5:][:]"../OUTPUT/T2.0/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 2.0",f20(x)

fit [20:] f10(x) "../OUTPUT/T1.0/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D10,a10
plot [20:][:]"../OUTPUT/T1.0/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 1.0",f10(x)

fit [30:] f08(x) "../OUTPUT/T0.8/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D08,a08
plot [30:][:]"../OUTPUT/T0.8/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.8",f08(x)

fit [80:] f07(x) "../OUTPUT/T0.7/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D07,a07
plot [50:][:]"../OUTPUT/T0.7/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.7",f07(x)

fit [50:] f06(x) "../OUTPUT/T0.6/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D06,a06
plot [50:][:]"../OUTPUT/T0.6/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.6",f06(x)

fit [200:] f055(x) "../OUTPUT/T0.55/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D055,a055
plot [200:][:]"../OUTPUT/T0.55/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.55",f055(x)

fit [250:] f052(x) "../OUTPUT/T0.52/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D052,a052
plot [250:][:]"../OUTPUT/T0.52/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.52",f052(x)

fit [1000:] f049(x) "../OUTPUT/T0.49/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D049,a049
plot [1000:][:]"../OUTPUT/T0.49/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.49",f049(x)

fit [7000:] f047(x) "../OUTPUT/T0.47/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D047,a047
plot [1000:][:]"../OUTPUT/T0.47/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.47",f047(x)

fit [10000:] f046(x) "../OUTPUT/T0.46/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D046,a046
plot [10000:][:]"../OUTPUT/T0.46/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.46",f046(x)

fit [10000:] f045(x) "../OUTPUT/T0.45/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) via D045,a045
plot [10000:][:]"../OUTPUT/T0.45/N1080/msd_NVT.txt" u 1:($2/(6*$1)):($3/(6*$1)) w errorl t "T = 0.45",f045(x)

Dfile="../THERMALIZE/data/D.txt"
system(sprintf("rm %s",Dfile))
set print Dfile
print "T D errD"
print 5.0,D50,D50_err
print 2.0,D20,D20_err
print 1.0,D10,D10_err
print 0.8,D08,D08_err
print 0.7,D07,D07_err
print 0.6,D06,D06_err
print 0.55,D055,D055_err
print 0.52,D052,D052_err
print 0.49,D049,D049_err
print 0.47,D047,D047_err
print 0.46,D046,D046_err
print 0.45,D045,D045_err
set print

reset
set out "./FIGURES/Dmsd.eps"
unset logs
set logs y
set ylabel "{/Times-Italic D}"
set xlabel "{/Times-Italic T}"
set ytics format "10^{%T}"
p Dfile u 1:2:3 w errorl not

set out "./FIGURES/Dmsd_Tc.eps"
reset
set logs
set key top left
set tics format "10^{%T}"
set ylabel "{/Times-Italic D}"
set xlabel "{/Times-Italic T-T}_c"
f(x)=a*(x-Tc)**expo
Tc=0.435; expo=1.75; a=.05
fit [:.7] f(x) Dfile u 1:2:3 via a,Tc,expo
p Dfile u ($1-Tc):($2):3 with errorl notitle, f(x+Tc) t "Power law fit"

