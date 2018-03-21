#!/bin/bash
#
# Script that runs ReadRun.py for all the samples.
# There is the possibility of running several samples in the same subscript,
# so that a same GPU can be loaded with more than one (useful with small systems).

#
# QUEUE PROPERTIES 
#
set -a #This must stay. It exports all the variables to ReadRun.sh (and other scripts)

readonly PROC_TAG="ra"
readonly USERNAME=`whoami`

readonly SYSTEM="PennPuter"
#readonly SYSTEM="talapas"; queue=gpu; simTime="0-02:00:00" #Days-HH:MM:SS

#PARAMETERS THAT SHOULD BE AT THE BEGINNING
hottestT=10.0
TLIST="2.0" # 0.6 0.49 0.466 0.44 0.43" #0.42
samLIST="0 1" # 1 2 3 4 5 6 7 8 9"
nsamples=`echo $samLIST|wc|awk '{print $2}'`

#This is to launch several runs on the same GPU
if [ $SYSTEM == "talapas" ]; then samBatchSize=2;
else samBatchSize=3; fi
numSamBatches=$(echo "($nsamples+$samBatchSize-1)/$samBatchSize"| bc) #workaround for ceil(nsamples/samBatchSize)
let numSamBatchesm1=$numSamBatches-1

#DIRECTORIES
scriptDIR=$PWD
rootDIR=$scriptDIR/../..
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES
mkdir -p $workDIR


#Range of temperatures
for T in $(echo $TLIST)
do
	cd $workDIR
	echo "T = $T"
	for Natoms in 65
	do
		echo "Natoms = $Natoms"
		inistateDIR=$workDIR/INITIAL-STATES/N$Natoms
		hottestTDIR=$workDIR/T$hottestT/N$Natoms
		paramsDIR=$workDIR/T$T/N$Natoms
		paramsFILE="$paramsDIR/params.in"

		mkdir -p $paramsDIR
		if ! [ -f "$paramsFILE" ]; then echo "$paramsFILE does not exist. Exiting..."; exit; fi

		for ibatch in $(seq 0 $numSamBatchesm1)
		do
			echo ibatch: $ibatch
			first=$(echo $ibatch*$samBatchSize|bc)
			last=$(echo "($ibatch+1)*$samBatchSize-1"|bc)
			fields=$(seq $first $last| awk '{a[NR]=$1}END{for(i=1;i<=NR;i++){printf "%d,",a[i]+1;}}' | sed 's/,$//')
			batch="$(echo $samLIST|cut -f$fields -d' ')"
			echo "batch:" $batch

			if [ $SYSTEM == 'talapas' ]; then
				sbatch $scriptDIR/ReadRun.sbatch
			elif [ $SYSTEM == 'PennPuter' ]; then
				# export batch
				bash $scriptDIR/ReadRun.sh #batch="$batch"

			else
				echo "SYSTEM = $SYSTEM was not recognised. Exiting..."
				exit

			fi



		done #ibatch

	done #N
done #T
