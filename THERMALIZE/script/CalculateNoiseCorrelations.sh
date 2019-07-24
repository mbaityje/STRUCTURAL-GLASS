#!/bin/bash

#SYSTEM
if [ `hostname` == "PennPuter" ] || [ `hostname` == "banshee" ];
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
LISTAT=${1:-"5.0"}
LISTAN=${2:-"1080"}
LISTATHERMOSTAT=${3:-"NVT"}

echo "lista T         : $LISTAT"
echo "lista N         : $LISTAN"
echo "lista thermostat: $LISTATHERMOSTAT"

shiftCFP="--shiftCFP"
showplots=""
showplots="--showplots"
kind='combined'
if [ -z $pot_mode ]; then pot_mode='xplor'; fi
if [ $maxtime   ]; then   maxtime="--maxtime=$maxtime"; fi
if [ $showplots ]; then showplots="--showplots"; fi
M=${M:-3}

declare -A TMIN=(  ["5.0"]=0.027 ["2.0"]=0.042 ["1.0"]=0.043 ["0.8"]=0.062 ["0.7"]=0.059 ["0.6"]=0.071 ["0.55"]=0.065 ["0.52"]=0.066 ["0.49"]=0.065 ["0.47"]=0.074 ["0.46"]=0.069 ["0.45"]=0.078 )
declare -A NCOEF=( ["5.0"]=9     ["2.0"]=7     ["1.0"]=7     ["0.8"]=5     ["0.7"]=5     ["0.6"]=7     ["0.55"]=7     ["0.52"]=6     ["0.49"]=6     ["0.47"]=5     ["0.46"]=6     ["0.45"]=7     )

for T in $(echo $LISTAT)
do
	for N in $(echo $LISTAN)
	do
		ncoef=${NCOEF[$T]}
		tmin=${TMIN[$T]}

		for thermostat in $(echo $LISTATHERMOSTAT)
		do
			cd $workDIR/T$T/N$N/
			echo "We are in $PWD"
			L=$(python $utilDIR/FindL.py ./S0/thermalized.gsd)
			if [ $noJK ]; 
			then
				python $exeDIR/CalculateNoiseCorrelationsLaplace.py -T$T --thermostat=$thermostat --kind=$kind --ncoef=$ncoef --M=$M --tmin=$tmin $showplots $shiftCFP $maxtime 
			else
				python $exeDIR/CalculateNoiseCorrelationsJK.py -T$T --thermostat=$thermostat -M=$M $maxtime --tmin=$tmin  --kind=$kind --ncoef=$ncoef --central_value='meanJK'
			fi
			echo ""
		done
	done
done

