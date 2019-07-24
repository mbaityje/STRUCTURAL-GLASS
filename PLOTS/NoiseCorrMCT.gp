#!/usr/bin/env gnuplot
#
# Plot the K(t) obtained through MCT
#                                                          #
############################################################
reset

kmax="_kmax28"
#kmax=""


set term post enh c eps
set out "FIGURES/K_MCT".kmax.".eps"

set logs x
set xlabel "{/Times-Italic t}"
set ylabel "{/Times-Italic K}({/Times-Italic t})"
set key top right samplen 1
set xtics format "10^{%T}"
load "scl.pal"
p "~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp ls 1 t "{/Times-Italic T} = 5.0  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T2.0/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp ls 2 t "{/Times-Italic T} = 2.0  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T1.0/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp t "{/Times-Italic T} = 1.0  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.8/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp t "{/Times-Italic T} = 0.8  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.7/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp t "{/Times-Italic T} = 0.7  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.6/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp t "{/Times-Italic T} = 0.6  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.55/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp t "{/Times-Italic T} = 0.55",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.52/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp t "{/Times-Italic T} = 0.52",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.49/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp t "{/Times-Italic T} = 0.49",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.47/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp t "{/Times-Italic T} = 0.47",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp t "{/Times-Italic T} = 0.46",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.45/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):3 w lp t "{/Times-Italic T} = 0.45"

set logs y
p "~/STRUCTURAL-GLASS/OUTPUT/T5.0/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp ls 1 t "{/Times-Italic T} = 5.0  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T2.0/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp ls 2 t "{/Times-Italic T} = 2.0  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T1.0/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp t "{/Times-Italic T} = 1.0  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.8/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp t "{/Times-Italic T} = 0.8  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.7/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp t "{/Times-Italic T} = 0.7  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.6/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp t "{/Times-Italic T} = 0.6  ",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.55/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp t "{/Times-Italic T} = 0.55",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.52/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp t "{/Times-Italic T} = 0.52",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.49/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp t "{/Times-Italic T} = 0.49",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.47/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp t "{/Times-Italic T} = 0.47",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp t "{/Times-Italic T} = 0.46",\
"~/STRUCTURAL-GLASS/OUTPUT/T0.45/N1080/KvertexFuchs_NVT".kmax.".txt" u ($1):2 w lp t "{/Times-Italic T} = 0.45"
set output
unset logs y


set term qt


reset
set logs x
set fit errorvariables

f045(x)=A045*exp(-(x/tau045)**beta045); beta045=0.6; tau045=200; A045=0.02
fit [1.5:14000] f045(x) "~/STRUCTURAL-GLASS/OUTPUT/T0.45/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 via A045,tau045,beta045
p [:]"~/STRUCTURAL-GLASS/OUTPUT/T0.45/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 w lp, f045(x)

f046(x)=A046*exp(-(x/tau046)**beta046); beta046=0.6; tau046=100; A046=0.01
fit [10:10000] f046(x) "~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 via A046,tau046,beta046
p [:]"~/STRUCTURAL-GLASS/OUTPUT/T0.46/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 w lp, f046(x)

f047(x)=A047*exp(-(x/tau047)**beta047); beta047=0.6; tau047=100; A047=0.01
fit [.6:2000] f047(x) "~/STRUCTURAL-GLASS/OUTPUT/T0.47/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 via A047,tau047,beta047
p [:]"~/STRUCTURAL-GLASS/OUTPUT/T0.47/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 w lp, f047(x)

f049(x)=A049*exp(-(x/tau049)**beta049); beta049=0.6; tau049=100; A049=0.01
fit [1.3:1700] f049(x) "~/STRUCTURAL-GLASS/OUTPUT/T0.49/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 via A049,tau049,beta049
p [:]"~/STRUCTURAL-GLASS/OUTPUT/T0.49/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 w lp, f049(x)

f052(x)=A052*exp(-(x/tau052)**beta052); beta052=0.6; tau052=100; A052=0.01
fit [1.:800] f052(x) "~/STRUCTURAL-GLASS/OUTPUT/T0.52/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 via A052,tau052,beta052
p [:]"~/STRUCTURAL-GLASS/OUTPUT/T0.52/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 w lp, f052(x)

f055(x)=A055*exp(-(x/tau055)**beta055); beta055=0.6; tau055=100; A055=0.01
fit [1.:300] f055(x) "~/STRUCTURAL-GLASS/OUTPUT/T0.55/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 via A055,tau055,beta055
p [:]"~/STRUCTURAL-GLASS/OUTPUT/T0.55/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 w lp, f055(x)

f06(x)=A06*exp(-(x/tau06)**beta06); beta06=0.6; tau06=10; A06=0.1
fit [2.:70] f06(x) "~/STRUCTURAL-GLASS/OUTPUT/T0.6/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 via A06,tau06,beta06
p [:]"~/STRUCTURAL-GLASS/OUTPUT/T0.6/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 w lp, f06(x)

f07(x)=A07*exp(-(x/tau07)**beta07); beta07=0.6; tau07=10; A07=0.05
fit [1.7:50] f07(x) "~/STRUCTURAL-GLASS/OUTPUT/T0.7/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 via A07,tau07,beta07
p [:]"~/STRUCTURAL-GLASS/OUTPUT/T0.7/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 w lp, f07(x)

f08(x)=A08*exp(-(x/tau08)**beta08); beta08=0.6; tau08=10; A08=0.1
fit [1.:18] f08(x) "~/STRUCTURAL-GLASS/OUTPUT/T0.8/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 via A08,tau08,beta08
p [:]"~/STRUCTURAL-GLASS/OUTPUT/T0.8/N1080/KvertexFuchs_NVT".kmax.".txt" u 1:3 w lp, f08(x)



tauFILE="../THERMALIZE/data/tauMCT".kmax.".txt"
set print tauFILE
print "T tau errtau beta errbeta A errA"
print  0.8, tau08, tau08_err, beta08, beta08_err, A08, A08_err
print  0.7, tau07, tau07_err, beta07, beta07_err, A07, A07_err
print  0.6, tau06, tau06_err, beta06, beta06_err, A06, A06_err
print 0.55,tau055,tau055_err,beta055,beta055_err,A055,A055_err
print 0.52,tau052,tau052_err,beta052,beta052_err,A052,A052_err
print 0.49,tau049,tau049_err,beta049,beta049_err,A049,A049_err
print 0.47,tau047,tau047_err,beta047,beta047_err,A047,A047_err
print 0.46,tau046,tau046_err,beta046,beta046_err,A046,A046_err
print 0.45,tau045,tau045_err,beta045,beta045_err,A045,A045_err
set print


# Plot tau
set logs
set xtics format "%g"
t(x)=a*(x-tc)**expo; expo=-2; tc=0.43; a=1
fit [:] t(x) tauFILE u 1:2:3 via a,tc,expo
p tauFILE u ($1-tc):2:3 w errorl, t(x+tc)


# Plot beta
unset logs
set xtics format "%g"
set xlabel "{/Times-Italic T}"
set ylabel "{/Symbol-Oblique b}"
p tauFILE u 1:4:5 w errorl

# Plot plateau
unset logs
set xtics format "%g"
set xlabel "{/Times-Italic T}"
set ylabel "Plateau"
p tauFILE u 1:6:7 w errorl


