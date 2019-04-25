#!/bin/bash


#SYSTEM
if [ `hostname` == "PennPuter" ] || [ `hostname` == "banshee" ] || [ `hostname` == "tango" ];
then readonly SCHEDULER="NONE";
elif [ -n $(hostname|grep talapas) ] || [ -n $(hostname|grep n[0-9][0-9][0-9]) ] || [ -n $(hostname|grep sn[0-9]) ]; 
then 
	readonly SCHEDULER='SLURM'
	queue=gpu
	simTime="1-00:00:00" #Days-HH:MM:SS
fi
readonly PROC_TAG="rt"
readonly USERNAME=$(whoami)


#DIRECTORIES
if [ $SCHEDULER == "SLURM" ];
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

ntw=${ntw:-1}
ntCd=${ntCd:-11}


echo "lista T         : $LISTAT"
echo "lista N         : $LISTAN"
echo "lista thermostat: $LISTATHERMOSTAT"




for T in $(echo $LISTAT)
do
	for N in $(echo $LISTAN)
	do
		for thermostat in $(echo $LISTATHERMOSTAT)
		do
			for ivecn in $(seq 1 $(wc -l $thermDIR/data/n.txt| cut -f1 -d" "))
			do
				n=$(sed "${ivecn}q;d" $thermDIR/data/n.txt)
				echo "n = $n"

				cd $workDIR/T$T/N$N/
				L="$(python $utilDIR/FindL.py ./S0/thermalized.gsd)"
				mkdir -p vertex
				cd vertex
				echo python $exeDIR/CalculateVertexCorrelations.py -L$L -T$T -N$N --thermostat=$thermostat -n $n --limit_input=$ntw --ntCd=$ntCd
				python $exeDIR/CalculateVertexCorrelations.py -L$L -T$T -N$N --thermostat=$thermostat -n $n --limit_input=$ntw --ntCd=$ntCd
			done
		done
	done
done

