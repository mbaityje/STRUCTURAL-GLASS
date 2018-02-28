#!/bin/bash


readonly PROC_TAG="cf"
readonly USERNAME=`whoami`

readonly SYSTEM="PennPuter"
#readonly SYSTEM="Talapas"; queue=gpu

#DIRECTORIES
scriptDIR=$PWD
thermDIR=$scriptDIR/..
rootDIR=$thermDIR/..
exeDIR=$thermDIR/progs
workDIR=$rootDIR/OUTPUT
mkdir -p $workDIR

#PARAMETERS THAT SHOULD BE AT THE BEGINNING
dt=0.0025
TLIST="0.466" # 0.6"
npoints=25
trajFreq=-$npoints
thermostat='NVT'

for T in $(echo $TLIST)
do
    if [ $T == 0.6 ]; then nsteps=`echo 400/$dt|bc`;
    elif [ $T == 0.466 ]; then nsteps=`echo 3000/$dt|bc`;
    fi
    
    cd $workDIR
    for Natoms in 1080
    do
	SAMPLEDIRS=`ls T$T/N$Natoms/S*/thermalized.gsd |sed 's:.*\/\(S.*\)\/.*:\1:'`

	for sam in $SAMPLEDIRS
	do
	    cd T$T/N$Natoms/$sam
	    if [ -f thermalized_backup.gsd ]; then conf=thermalized_backup.gsd
	    else conf=thermalized.gsd;
	    fi
	    echo "sample: $sam"
	    python $exeDIR/corrF.py --user="$conf -N$Natoms --addsteps=True -t$nsteps --trajFreq=$trajFreq --thermostat=$thermostat --dt=$dt --temperature=$T"
	    cd $workDIR
	done
    done
done


