#!/bin/bash

#Process tag is used for the queue scheduler
readonly PROC_TAG="cs"

#Script Options
doTS=0 #0: only does IS, 1: does IS and TS

readonly SYSTEM="PennPuter"
#readonly SYSTEM="Talapas"

#DIRECTORIES
thermDIR=$PWD/..
rootDIR=$PWD/../..
exeDIR=$thermDIR/progs
scriptDIR=$thermDIR/script
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES


#Parameters
readonly kB=0.00831445986144858 #Boltzman constant in our units
readonly dt=0.0025
readonly ttot=`echo 10^7|bc` #We will want 10^9 steps
TLIST="0.43" #"10.0 2.0 0.6 0.466 0.44 0.43 0.42 0.41"

cd $workDIR
for T in $(echo $TLIST)
do
	echo "T = $T"
	cd T$T
	deltaE=0.5
#	deltaE=0.2
#	deltaE=$(echo 10*$kB*$T | bc -lq)

	for Natoms in 65
	do
		cd N$Natoms
		for ISAM in $(ls -d S* | sed 's/S//') #ISAM=any sample that has been simulated
		do
			cd S$ISAM
			echo ISAM = $ISAM
			pwd
			bash $scriptDIR/ChunkSect.sh $T $dt $deltaE $ttot $doTS
			cd ..
		done #ISAM
		cd ..
	done #Natoms
	cd ..
done #T
