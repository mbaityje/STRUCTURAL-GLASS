#!/bin/bash

#Process tag is used for the queue scheduler
readonly PROC_TAG="cs"

#Script Options
readonly doridge=1 #0: only does IS, 1: does IS and Ridge

if [ `hostname` == "PennPuter" ] || [ `hostname` == "tango" ];
then
	SYSTEM="PennPuter";
else
	SYSTEM="talapas";
	queue=gpu
	simTime="0-12:00:00" #Days-HH:MM:SS
fi

#DIRECTORIES
rootDIR=$PWD/../..
thermDIR=$rootDIR/THERMALIZE
scriptDIR=$thermDIR/script
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES


#Parameters
readonly dt=0.0025
readonly ttot=`echo 10^9|bc` #We will want 10^9 steps
TLIST="0.49" #"10.0 2.0 0.6 0.466 0.44 0.43 0.42 0.41"
readonly deltaE=0.0001

cd $workDIR
for T in $(echo $TLIST)
do
	echo "T = $T"
	cd T$T


	for Natoms in 65
	do
		cd N$Natoms
		for ISAM in 0 #$(ls -d S* | sed 's/S//') #ISAM=any sample that has been simulated
		do
			cd S$ISAM
			echo ISAM = $ISAM
			pwd
			if [ $SYSTEM == "PennPuter" ]
			then
				bash $scriptDIR/ChunkSect.sh $T $dt $deltaE $ttot $doridge
			elif [ $SYSTEM == "talapas" ]
			then
				nombre=N$Natoms${PROC_TAG}T${T}i${ISAM}
				if [ 0 == `squeue -u$(whoami) -n $nombre|grep $(whoami)|wc -l` ]
				then
				echo "sbatch --job-name=$nombre --export=T=$T,dt=$dt,deltaE=$deltaE,ttot=$ttot,doridge=$doridge $scriptDIR/ChunkSect.sh"
				sbatch  -p $queue --time=$simTime --job-name=$nombre --export=T=$T,dt=$dt,deltaE=$deltaE,ttot=$ttot,doridge=$doridge $scriptDIR/ChunkSect.sh
				fi
			fi
			cd ..
		done #ISAM
		cd ..
	done #Natoms
	cd ..
done #T
