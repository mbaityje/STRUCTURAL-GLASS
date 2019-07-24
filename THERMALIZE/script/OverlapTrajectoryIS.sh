#!/bin/bash


rootDIR=$PWD/../..
thermDIR=$rootDIR/THERMALIZE
scriptDIR=$thermDIR/script
dataDIR=$thermDIR/data
progDIR=$thermDIR/progs
outDIR=$rootDIR/OUTPUT

ntbar=100

for T in 5.0 2.0 1.0 0.8 0.7 0.6 0.55 0.52 0.49
do
    for N in 65
    do
	for ISAM in 0 1 2 3 4 5 6 7 8 9
	do
	    samDIR=$outDIR/T$T/N$N/shift/S$ISAM
	    cd $samDIR
	    inConf=thermalized.gsd
	    for irep in {0..9}
	    do
			nsteps=$(awk -vT=$T -vN=$N '($1==T && $2==N){printf("%d\n",$3/0.0025)}' $dataDIR/tau_alpha.txt)

			echo "nsteps = $nsteps"
			seed=$RANDOM
			python $progDIR/OverlapTrajectoryIS.py --user="$inConf --nsteps=$nsteps --temperature=$T --ntbar=$ntbar --irep=$irep --seed=$seed --iframe=0" &
	    done #irep
	    wait
	done #isam
    done #N
done #T
