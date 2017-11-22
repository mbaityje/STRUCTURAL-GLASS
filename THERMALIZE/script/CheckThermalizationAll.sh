#!/bin/bash
#
# Check thermalization of all the samples by calculating
# the self-intermediate scattering function, and tau.
#

readonly PROC_TAG="rt"

#Script Options
tau_of_t=0 #1: calculate Fkt on all the heavyTraj, 0: calculate Fkt on only the last configuration


readonly SYSTEM="PennPuter"
#readonly SYSTEM="Talapas"

#DIRECTORIES
scriptDIR=$PWD
thermDIR=$PWD/..
rootDIR=$PWD/../..
exeDIR=$thermDIR/progs
scriptDIR=$thermDIR/script
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES
mkdir -p $workDIR
cd $workDIR



#Each T requires a different nsteps
readonly dt=0.0025
declare -A NSTEPS_LIST=( ["10.0"]=$(echo 0.1/$dt |bc) ["2.0"]=$(echo 1.0/$dt |bc) )




for Tdir in `ls -d T*`
do
    T=`echo $Tdir | sed 's/^T//'`
	echo "T=$T"
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
	    	if [ 1 -eq $tau_of_t ]
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
    			bash $scriptDIR/SelfIntermediateScatteringFunction.sh $thermalizeFile 0 $nsteps $T $dt $tau_of_t
    		fi
	    	cd ..
		done
		cd ..
    done
    cd ..
done

