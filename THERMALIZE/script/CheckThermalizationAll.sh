#!/bin/bash
# Lines for slurm
#SBATCH --ntasks=1
#SBATCH -p longgpu # partition (queue) 
#SBATCH --gres=gpu:1  
#
# Check thermalization of all the samples by calculating
# the self-intermediate scattering function, and tau.
#

readonly PROC_TAG="ct"
readonly USERNAME=`whoami`

#Script Options
tau_of_t=0 #1: calculate Fkt on all the heavyTraj, 0: calculate Fkt on only the last configuration


#readonly SYSTEM="PennPuter"
readonly SYSTEM="Talapas"

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
declare -A NSTEPS_LIST=( ["10.0"]=$(echo 10/$dt |bc) ["2.0"]=$(echo 100.0/$dt |bc) ["0.6"]=$(echo 2000.0/$dt |bc)  ["0.466"]=$(echo 20000.0/$dt |bc) ["0.44"]=$(echo 40000.0/$dt |bc) ["0.43"]=$(echo 80000.0/$dt |bc) ["0.42"]=$(echo 160000.0/$dt |bc) ["0.41"]=$(echo 400000.0/$dt |bc))




for Tdir in T2.0 #`ls -d T*|sort -r`
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
	cd $Ndir
	for SAMdir in `ls -d S*`
	do
	    ISAM=`echo $SAMdir | sed 's/^S//'`
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
		elif [ $SYSTEM == "Talapas" ]
		then
		    nombre=N$Natoms${PROC_TAG}T${T}i${isam}
		    if [ 0 == `squeue -u$USERNAME -n $nombre|grep $USERNAME|wc -l` ]
		    then
    			sbatch --job-name=$nombre --export=filename=$thermalizeFile,iframe=0,nsteps=$nsteps,T=$T,dt=$dt,tau_of_t=$tau_of_t $scriptDIR/SelfIntermediateScatteringFunction.sh
		    fi
		fi
    	    fi
	    cd ..
	done #ISAM
	cd ..
    done #N
    cd ..
done #T

