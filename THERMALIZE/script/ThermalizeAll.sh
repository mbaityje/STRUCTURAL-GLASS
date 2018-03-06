#!/bin/bash
#
#A partire dalle configurazioni calde ne creo di piu` fredde

#
# QUEUE PROPERTIES 
#
readonly PROC_TAG="rt"
readonly USERNAME=`whoami`

#readonly SYSTEM="PennPuter"
readonly SYSTEM="talapas"; queue=gpu; simTime="0-02:00:00" #Days-HH:MM:SS

#PARAMETERS THAT SHOULD BE AT THE BEGINNING
nsam=10
let nsamm1=$nsam-1
dt=0.0025
backupFreq=`echo 10/$dt|bc`
hottestT=10.0
TLIST="0.6" #"2.0 0.6 0.49 0.466 0.44 0.43" #0.42

#DIRECTORIES
scriptDIR=$PWD
thermDIR=$PWD/..
rootDIR=$PWD/../..
exeDIR=$thermDIR/progs
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES
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
	    echo "isam: $isam"
	    mkdir -p S$isam
	    cd S$isam
	    
	    seed=$(od -vAn -N4 -tu4 < /dev/urandom)
	    thermConfName=thermalized.gsd
	    case $T in
		10.0)  totMDsteps=$(echo 0.5*10^3/$dt|bc); thermostat=MB;  tauT=1.0; heavyTrajFreq=`echo $totMDsteps/4|bc`;   queue=short;   initConf=$inistateDIR/initIS.gsd;;
		2.0)   totMDsteps=$(echo 1*10^4/$dt  |bc); thermostat=NVT; tauT=0.1; heavyTrajFreq=`echo $totMDsteps/4|bc`;   queue=short;   initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.6)   totMDsteps=$(echo 3*10^6/$dt  |bc); thermostat=NVT; tauT=0.1; heavyTrajFreq=`echo $totMDsteps/20|bc`;  queue=short;   initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.49)  totMDsteps=$(echo 1*10^8/$dt  |bc); thermostat=NVT; tauT=0.1; heavyTrajFreq=`echo $totMDsteps/40|bc`;  queue=standard;initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.466) totMDsteps=$(echo 1*10^9/$dt  |bc); thermostat=NVT; tauT=0.1; heavyTrajFreq=`echo $totMDsteps/50|bc`;  queue=gpu;     initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.44)  totMDsteps=$(echo 5*10^9/$dt  |bc); thermostat=NVT; tauT=0.1; heavyTrajFreq=`echo $totMDsteps/100|bc`; queue=gpu;     initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.43)  totMDsteps=$(echo 5*10^10/$dt |bc); thermostat=NVT; tauT=0.1; heavyTrajFreq=`echo $totMDsteps/100|bc`; queue=gpu;     initConf=$hottestTDIR/S$isam/$thermConfName;;
		0.42)  totMDsteps=$(echo 1*10^11/$dt |bc); thermostat=NVT; tauT=0.1; heavyTrajFreq=`echo $totMDsteps/100|bc`; queue=gpu;     initConf=$hottestTDIR/S$isam/$thermConfName;;
		*) echo "How many steps for this temperature?";exit;;
	    esac
	    
	    #If thermalized configuration already exists, then continue from there (to be able to make longer thermalization runs)
	    if [ -e $thermConfName ] ;then
		initConf=$thermConfName
		initConfSteps=$(python $utilDIR/FindNsteps.py $initConf)
		echo "Number of steps in $initConf: $initConfSteps"
		echo "Target number of steps: $totMDsteps"
		if [ $initConfSteps -ge $totMDsteps ]; then
		    echo "No more steps are needed"
		    cd ..; continue
		fi
	    else
		echo "Starting run from $initConf"
	    fi
	    
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
		echo `whoami`$USERNAME@`uname -n` > thermalized.time
		time (python $exeDIR/ReadAndThermalize.py --user="$initConf -N$Natoms -s$seed -T$T -t$totMDsteps --tauT=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=$backupFreq --heavyTrajFreq=$heavyTrajFreq"  2>&1) 2>>thermalized.time
	    elif [ $SYSTEM == "talapas" ]; 
	    then
		nombre=N${Natoms}${PROC_TAG}T${T}i${isam}
		if [ 0 == `squeue -u$USERNAME -n $nombre|grep $USERNAME|wc -l` ]
		then
		    echo sbatch --job-name=$nombre -p $queue --export=exeDIR=$exeDIR,initConf=$initConf,Natoms=$Natoms,seed="${seed}",T=$T,totMDsteps=$totMDsteps,tauT=$tauT,dt=$dt,thermostat=$thermostat,backupFreq=$backupFreq,heavyTrajFreq=$heavyTrajFreq $scriptDIR/Thermalize.sbatch
		    sbatch --job-name=$nombre -p $queue --time=$simTime --export=exeDIR=$exeDIR,initConf=$initConf,Natoms=$Natoms,seed="${seed}",T=$T,totMDsteps=$totMDsteps,tauT=$tauT,dt=$dt,thermostat=$thermostat,backupFreq=$backupFreq,heavyTrajFreq=$heavyTrajFreq $scriptDIR/Thermalize.sbatch
		fi
	    else
		echo "SYSTEM=$SYSTEM not recognized"
		exit
	    fi
	    echo ""
	    cd ..
	done #isam
	cd ..
    done #Natoms
    cd ..
    wait
done #T
cd ..
