#!/bin/bash
#A partire dalle configurazioni calde ne creo di piu` fredde

#
# QUEUE PROPERTIES 
#
readonly PROC_TAG="rt"

readonly SYSTEM="PennPuter"
#readonly SYSTEM="Talapas"
QUEUE="gpu"
walltime="0:23:30:00"

#PARAMETERS THAT SHOULD BE AT THE BEGINNING
nsam=1
let nsamm1=$nsam-1
dt=0.0025
backupFreq=`echo 10/$dt|bc`
hottestT=10.0
TLIST="10.0 2.0" #"10.0 2.0 0.6 0.466 0.44 0.43 0.42 0.41"

#DIRECTORIES
scriptDIR=$PWD
thermDIR=$PWD/..
rootDIR=$PWD/../..
exeDIR=$thermDIR/progs
workDIR=$rootDIR/OUTPUT
mkdir -p $workDIR
cd $workDIR




#Range of temperatures
for T in $(echo $TLIST)
do
    mkdir -p T$T
    cd T$T
    echo "T = $T"

    for Natoms in 65
    do
	mkdir -p N$Natoms
	cd N$Natoms
	echo "Natoms = $Natoms"

	inistateDIR=$workDIR/INITIAL-STATES/N$Natoms
	hottestTDIR=$workDIR/T$hottestT/N$Natoms

	
	for isam in $(seq 0 $nsamm1)
	do
	    mkdir -p S$isam
	    cd S$isam
		
	    seed=$(od -vAn -N4 -tu4 < /dev/urandom) #seed is stored in measurement file
	    thermConfName=thermalized.gsd
	    case $T in
		10.0)  totMDsteps=8000000;     thermostat=MB;  tau=1.0; heavyTrajFreq=`echo $totMDsteps/4|bc`; initConf=$inistateDIR/initIS.gsd;;
		2.0)   totMDsteps=10000000;    thermostat=NVT; tau=0.1; heavyTrajFreq=`echo $totMDsteps/4|bc`; initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.6)   totMDsteps=20000000;    thermostat=NVT; tau=0.1; heavyTrajFreq=`echo $totMDsteps/4|bc`; initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.466) totMDsteps=40000000;    thermostat=NVT; tau=0.1; heavyTrajFreq=`echo $totMDsteps/4|bc`; initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.44)  totMDsteps=400000000;   thermostat=NVT; tau=0.1; heavyTrajFreq=`echo $totMDsteps/4|bc`; initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.43)  totMDsteps=2000000000;  thermostat=NVT; tau=0.1; heavyTrajFreq=`echo $totMDsteps/4|bc`; initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.42)  totMDsteps=4000000000;  thermostat=NVT; tau=0.1; heavyTrajFreq=`echo $totMDsteps/4|bc`; initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.41)  totMDsteps=10000000000; thermostat=NVT; tau=0.1; heavyTrajFreq=`echo $totMDsteps/4|bc`; initConf=$hottestTDIR/S$isam/$thermConfName;;
		*) echo "How many steps for this temperature?";exit;;
	    esac
	    
	    echo $initConf
	    
	    #If thermalized configuration already exists, then continue from there (to be able to make longer thermalization runs)
	    if [ -e $thermConfName ] ;then
			initConf=$thermConfName
	    fi

	    echo INITCONF: $(ls $PWD/$initConf)
	    
	    #It may happen that no initial configuration is available.
	    #In that case, give a warning and skip sample
	    if [ ! -f $initConf ]; then
			echo "WARNING: $initconf does not exist"
			cd ..; continue
	    fi
	
	    #
	    # Run Python script
	    #
	    if [ $SYSTEM == "PennPuter" ]; then
			echo `whoami`@`uname -n` > thermalized.time
			time (python $exeDIR/ReadAndThermalize.py --user="$initConf -N$Natoms -s$seed -T$T -t$totMDsteps --tau=$tau --dt=$dt --thermostat=$thermostat --backupFreq=$backupFreq --heavyTrajFreq=$heavyTrajFreq"  2>&1) 2>>thermalized.time
	    else
			echo "SYSTEM=$SYSTEM not recognized"
		exit
	    fi
	    cd ..
	done #isam
	cd ..
    done #Natoms
    cd ..
    wait
done #T
cd ..
