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
if [ -z $dt       ]; then       dt=0.0025 ; fi
if [ $maxtime   ]; then   maxtime="--maxtime=$maxtime"; fi
if [ $softening ]; then softening='--softening'       ; fi
if [ $tstar     ]; then     tstar="--tstar=$tstar"    ; fi
if [ $normalsc  ]; then  normalsc='--normalsc'        ; fi
if [ $lin       ]; then       lin='--lin'             ; fi
if [ $linsc     ]; then     linsc='--linlsc'          ; fi
if [ $fits      ]; then      fits='--fits'            ; fi


for T in $(echo $LISTAT)
do
	for N in $(echo $LISTAN)
	do
		M=3
		if [ $T == 0.7 ];
		then
		    ncoef=10
		elif [ $T == 0.6 ];
		then
		     ncoef=9
		elif [ $T == 0.55 ];
		then
		     ncoef=12
		elif [ $T == 0.49 ];
		then
		     ncoef=9
		elif [ $T == 0.47 ];
		then
		     ncoef=7
		elif [ $T == 0.46 ];
		then
		     ncoef=8
		else
		    ncoef=10
		fi
		for thermostat in $(echo $LISTATHERMOSTAT)
		do
			cd $workDIR/T$T/N$N/
			echo "We are in $PWD"
			L=$(python $utilDIR/FindL.py ./S0/thermalized.gsd)
			# python $exeDIR/CalculateNoiseCorrelations.py -L$L -T$T -N$N --dt=$dt --thermostat=$thermostat $shiftCFP $maxtime $softening $tstar $normalsc $lin $linsc $fits
			echo "python $exeDIR/CalculateNoiseCorrelationsLaplace.py -L$L -T$T -N$N --dt=$dt --thermostat=$thermostat --kind=$kind --M=$M --showplots=$showplots $shiftCFP $maxtime $softening $tstar $normalsc $lin $linsc $fits"
			python $exeDIR/CalculateNoiseCorrelationsLaplace.py -L$L -T$T -N$N --dt=$dt --thermostat=$thermostat --kind=$kind --ncoef=$ncoef --M=$M $showplots $shiftCFP $maxtime $softening $tstar $normalsc $lin $linsc $fits
			echo ""
		done
	done
done

