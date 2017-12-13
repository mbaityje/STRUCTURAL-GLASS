#!/bin/bash


rootDIR=$PWD/../..
thermDIR=$rootDIR/THERMALIZE
scriptDIR=$thermDIR/script
progDIR=$thermDIR/progs
outDIR=$rootDIR/OUTPUT

declare -A NSTEPS_LIST=( ["10.0"]=1000 ["2.0"]=4000 ["0.6"]=0  ["0.466"]=0 ["0.44"]=0 ["0.43"]=$(echo 10^8|bc) ["0.42"]=0 ["0.41"]=0 )
ntbar=100

for T in 10.0 2.0 0.43
do
    for N in 65
    do
	for ISAM in 0
	do
	    samDIR=$outDIR/T$T/N$N/S$ISAM
	    cd $samDIR
	    inConf=thermalized.gsd
	    for irep in {0..4}
	    do
		nsteps=${NSTEPS_LIST[$T]}
		echo "nsteps = $nsteps"
		seed=$RANDOM
		python $progDIR/OverlapTrajectoryIS.py --user="$inConf --nsteps=$nsteps --temperature=$T --ntbar=$ntbar --irep=$irep --seed=$seed --iframe=0" &
		
		
	    done #irep
	    wait
	done #isam
    done #N
done #T
