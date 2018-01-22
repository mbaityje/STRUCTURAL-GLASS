#!/bin/bash
#
# Check thermalization of all the samples by calculating
# the self-intermediate scattering function, and tau.
#

readonly PROC_TAG="ct"
readonly USERNAME=`whoami`
queue=longgpu
#queue=short


#Script Options
readonly tau_of_t=0 #1: calculate Fkt on all the heavyTraj, 0: calculate Fkt on only the last configuration

if [ `hostname` == "PennPuter" ]; then SYSTEM="PennPuter";
else SYSTEM="talapas"; fi

#DIRECTORIES
rootDIR=$PWD/../..
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
scriptDIR=$thermDIR/script
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES
mkdir -p $workDIR
cd $workDIR

echo "Adesso mi trovo in $PWD"


#Each T requires a different nsteps
readonly dt=0.0025
#Of these nsteps, the following are fine tuned: T=10.0,2.0
#Questi sono per |k|=46
#declare -A NSTEPS_LIST=( ["10.0"]=$(echo 0.1/$dt |bc) ["2.0"]=$(echo 0.2/$dt |bc) ["0.6"]=$(echo 0.5/$dt |bc)  ["0.49"]=$(echo 2.0/$dt |bc) ["0.466"]=$(echo 8.0/$dt |bc) ["0.44"]=$(echo 16.0/$dt |bc) ["0.43"]=$(echo 32.0/$dt |bc) ["0.42"]=$(echo 64.0/$dt |bc) ["0.41"]=$(echo 128.0/$dt |bc))
declare -A NSTEPS_LIST=( ["10.0"]=$(echo 0.5/$dt |bc) ["2.0"]=$(echo 5.0/$dt |bc) ["0.6"]=$(echo 100.0/$dt |bc) ["0.43"]=$(echo 320.0/$dt |bc) )




for Tdir in `ls -d T*|sort -r`
do
    T=`echo $Tdir | sed 's/^T//'`

    echo "-------------------------------------------"
    echo "|xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|"
    echo "-------------------------------------------"
    echo "T=$T"
    pwd
    cd $Tdir
    nsteps=${NSTEPS_LIST[$T]}
    
    for Ndir in `ls -d N*`
    do
	N=`echo $Ndir | sed 's/^N//'`
	echo "N=$N"
	cd $Ndir
	for SAMdir in `ls -d S*`
	do
	    ISAM=`echo $SAMdir | sed 's/^S//'`
	    echo "ISAM=$ISAM"
	    cd $SAMdir
	    if [ 1 -eq $tau_of_t ] #I never tried this case
    	    then
		heavyTrajFile=heavyTraj.gsd
		Nframes=`python $utilDIR/FindNFrames.py $heavyTrajFile`
	    	let Nframesm1=$Nframes-1
		for iframe in $(seq 0 $Nframesm1)
		do
		    bash $scriptDIR/SelfIntermediateScatteringFunction.sh $heavyTrajFile $iframe $nsteps $T $dt $tau_of_t
	    	done
    	    else
    		thermalizeFile=thermalized.gsd
		if ! [ -f $thermalizeFile ]; then echo "$PWD/$thermalizeFile does not exist"; cd ..; continue; fi
		
		if [ $SYSTEM == "PennPuter" ]
		then
    		    bash $scriptDIR/SelfIntermediateScatteringFunction.sh $thermalizeFile 0 $nsteps $T $dt $tau_of_t
		elif [ $SYSTEM == "talapas" ]
		then
		    nombre=N$N${PROC_TAG}T${T}i${ISAM}
		    if [ 0 == `squeue -u$USERNAME -n $nombre|grep $USERNAME|wc -l` ]
		    then
    			sbatch --job-name=$nombre -p$queue --export=filename=$thermalizeFile,iframe=0,nsteps=$nsteps,T=$T,dt=$dt,tau_of_t=$tau_of_t $scriptDIR/SelfIntermediateScatteringFunction.sh
		    fi
		fi
    	    fi
	    cd ..
	done #ISAM
	cd ..
    done #N
    cd ..
done #T

