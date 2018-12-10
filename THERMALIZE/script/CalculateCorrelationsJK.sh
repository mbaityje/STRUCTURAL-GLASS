#!/bin/bash


#SYSTEM
if [ `hostname` == "PennPuter" ] || [ `hostname` == "banshee" ];
then SYSTEM="PennPuter";
elif [ `hostname` == "tango" ];
then SYSTEM="PennPuter";
else SYSTEM="talapas"; fi
echo SYSTEM = $SYSTEM


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
cd $workDIR
thermtimesFILE=$dataDIR/thermalizationtimes.txt

#PARAMETERS
#Parameters are taken from command line, though I set some default value, mainly for code developing and debugging
observables=${1:-"--msd --Fkt --CFF --CFP --CPP"}
LISTAT=${2:-"5.0"}
LISTAN=${3:-"1080"}
LISTATHERMOSTAT=${4:-"NVT"}

if [ $limit_input ]
then limit_input="--limit_input=$limit_input"
else limit_input=""
fi


echo "lista T         : $LISTAT"
echo "lista N         : $LISTAN"
echo "lista thermostat: $LISTATHERMOSTAT"

readonly dt=0.0025
pot_mode='xplor'
trajFreq=-1000 #A negative trajFreq means we sample |trajFreq| logarithmically distributed bins



nblo=10 #Only 10 JK blocks because it takes long to calculate K(t) on each single block

for T in $(echo $LISTAT)
do
	for N in $(echo $LISTAN)
	do
		for thermostat in $(echo $LISTATHERMOSTAT)
		do
			cd $workDIR/T$T/N$N/
			L="$(python $utilDIR/FindL.py ./S0/thermalized.gsd)"
			ntrajtot=$(ls S*/trajectories/trajectory${thermostat}*.gsd|wc -l) 
			lblo=$(echo "$ntrajtot/$nblo"|bc) 
			echo "lblo = $lblo"
			echo "python $exeDIR/CalculateCorrelationsJK.py -L$L -T$T -N$N --dt=$dt --thermostat=$thermostat --lblo=$lblo $observables $limit_input"
			python $exeDIR/CalculateCorrelationsJK.py -L$L -T$T -N$N --dt=$dt --thermostat=$thermostat --lblo=$lblo $observables $limit_input
		done
	done
done

