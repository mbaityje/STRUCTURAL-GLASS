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

if [ -z $pot_mode ]; then pot_mode='xplor'; fi
if [ -z $dt       ]; then       dt=0.0025 ; fi
if [ $shiftCFP  ]; then  shiftCFP='--shiftCFP'        ; fi
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
		for thermostat in $(echo $LISTATHERMOSTAT)
		do
			cd $workDIR/T$T/N$N/
			ls $PWD
			L=$(python $utilDIR/FindL.py ./S0/thermalized.gsd)
			echo python $exeDIR/CalculateNoiseCorrelations.py -L$L -T$T -N$N --dt=$dt --thermostat=$thermostat $shiftCFP $maxtime $softening $tstar $fits
			python $exeDIR/CalculateNoiseCorrelations.py -L$L -T$T -N$N --dt=$dt --thermostat=$thermostat $shiftCFP $maxtime $softening $tstar $normalsc $lin $linsc $fits
		done
	done
done

