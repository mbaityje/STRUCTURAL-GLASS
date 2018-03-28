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

#readonly SYSTEM='PennPuter'
readonly SYSTEM='talapas'

if [ $SYSTEM == 'talapas' ]; then
    # enforce python 3
    module unload anaconda2
    module load anaconda3
    queue=gpu
    simTime="0-02:00:00" #Days-HH:MM:SS
fi



#PARAMETERS THAT SHOULD BE AT THE BEGINNING
hottestT=10.0
TLIST="0.6 0.49 0.466 0.44 0.43"
samLIST="0 1 2 3 4 5 6 7 8 9"
nsamples=`echo $samLIST|wc|awk '{print $2}'`



#This is to launch several runs on the same GPU (in the case we are on talapas and the queue is gpu)
if [ $SYSTEM == "talapas" ]; 
then 
    #If we run on gpus, make batches of samples to run several 
    #on the same gpu
    if [ `echo $queue|grep gpu` == '' ]; then 
	samBatchSize=1;
    else
	samBatchSize=3;
    fi    
else samBatchSize=1
fi
numSamBatches=$(echo "($nsamples+$samBatchSize-1)/$samBatchSize"| bc) #workaround for ceil(nsamples/samBatchSize)
let numSamBatchesm1=$numSamBatches-1

#DIRECTORIES
scriptDIR=$PWD
rootDIR=$scriptDIR/../..
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES
logDIR=$rootDIR/LOGS
mkdir -p $workDIR


#Range of temperatures
for T in $(echo $TLIST)
do
	cd $logDIR
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

			if [ $SYSTEM == 'talapas' ]; 
			then
			    nombre=N${Natoms}${PROC_TAG}T${T}b${ibatch}
			    if [ 0 == `squeue -u$USERNAME -n $nombre|grep $USERNAME|wc -l` ]
			    then
				sbatch --job-name=$nombre -p $queue --time=$simTime $scriptDIR/ReadRun.sh
			    fi
			elif [ $SYSTEM == 'PennPuter' ];
			then
			    bash $scriptDIR/ReadRun.sh
			else
				echo "SYSTEM = $SYSTEM was not recognised. Exiting..."
				exit
			fi
		done #ibatch
	done #N
done #T
