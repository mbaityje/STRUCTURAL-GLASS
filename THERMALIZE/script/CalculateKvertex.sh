#!/bin/bash


SYSTEM="PennPuter"

#DIRECTORIES
if [ $SYSTEM == "talapas" ];
then rootDIR=/home/mbaity/STRUCTURAL-GLASS/
else 
	rootDIR=$PWD/../..
	echo "Setting rootDIR: $rootDIR"
fi
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
dataDIR=$thermDIR/data
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES

LISTAT=${1:-"5.0 0.46"}
LISTAN=${2:-"1080"}
LISTATHERMOSTAT=${3:-"NVT"}
if [ $showplots ]; then showplots="--showplots"; fi
if [ $kmax      ]; then kmax="--kmax=$kmax"; fi
if [ $nofuchs   ]; then fuchs=""; else fuchs="--fuchs"; fi


for T in $(echo $LISTAT)
do
	echo "T = $T"
	for N in $(echo $LISTAN)
	do
		for thermostat in $(echo $LISTATHERMOSTAT)
		do	
			cd $workDIR/T$T/N$N/vertex/
			L=$(python $utilDIR/FindL.py ../S0/thermalized.gsd)
			pwd
			echo "python $exeDIR/CalculateKvertex.py -T$T -L$L -N$N --thermostat=$thermostat $kmax $showplots $fuchs"
			python $exeDIR/CalculateKvertex.py -T$T -L$L -N$N --thermostat=$thermostat $kmax $showplots $fuchs
		done
	done
done

