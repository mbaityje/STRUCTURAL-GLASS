#!/bin/bash
#
#A partire dalle configurazioni calde ne creo di piu` fredde

#
# QUEUE PROPERTIES 
#
readonly PROC_TAG="rt"
readonly USERNAME=`whoami`

readonly SYSTEM="PennPuter"
#readonly SYSTEM="talapas"

queue=gpu

#PARAMETERS THAT SHOULD BE AT THE BEGINNING
nsam=10
let nsamm1=$nsam-1
dt=0.0025
backupFreq=`echo 10/$dt|bc`
hottestT=10.0
TLIST=${1:-"5.0 2.0 1.0 0.8 0.6"}

#DIRECTORIES
scriptDIR=$PWD
thermDIR=$PWD/..
rootDIR=$PWD/../..
exeDIR=$thermDIR/progs
dataDIR=$thermDIR/data
workDIR=$rootDIR/OUTPUT
mkdir -p $workDIR
cd $workDIR
thermtimesFILE=$dataDIR/thermalizationtimes.txt


#Range of temperatures
for T in $(echo $TLIST)
do
	mkdir -p T$T
	cd T$T
	echo "T = $T"
	
	for Natoms in 1080
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

		#Set parameters. Hottest T is done with Andersen Thermostat, so it has different parameters
		thermConfName=thermalized.gsd
		if [ $T == $hottestT ];
		then
			totMDsteps=$(echo 10/$dt|bc)
			thermostat=MB
			heavyTrajFreq=0
			tauT=1.0
			initConf=$inistateDIR/initIS.gsd
		else
			tautherm=$(awk -vT=$T -vN=$Natoms '($1==T && $2==N){print $3}' $thermtimesFILE)
			tauthermsteps=$(echo $tautherm/$dt|bc)
			totMDsteps=$(echo 10*$tauthermsteps|bc)
			thermostat=NVT
			heavyTrajFreq=$tauthermsteps
			tauT=0.1
			initConf=$hottestTDIR/S$isam/$thermConfName
		fi

		#If (partially) thermalized configuration already exists, then continue from there (to be able to make longer thermalization runs)
		echo Proposed initConf: $initConf
		if [ -e $thermConfName ] ;then
			initConf=$thermConfName
			echo "Continuing run from $initConf"
			startfromzero=''
		else
			echo "Starting run from $initConf"
			startfromzero='--startfromzero'
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
		# BEWARE: there is a & at the end of the following command, which means that samples will be run in parallel
		time (python $exeDIR/ReadAndThermalize.py --user="$initConf -N$Natoms -s$seed -T$T -t$totMDsteps --tauT=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=$backupFreq --heavyTrajFreq=$heavyTrajFreq $startfromzero"  2>&1) 2>>thermalized.time #&
		elif [ $SYSTEM == "talapas" ]; 
		then
		nombre=N${Natoms}${PROC_TAG}T${T}i${isam}
		if [ 0 == `squeue -u$USERNAME -n $nombre|grep $USERNAME|wc -l` ]
		then
			echo "sbatch --job-name=$nombre --export=exeDIR=$exeDIR,initConf=$initConf,Natoms=$Natoms,seed="${seed}",T=$T,totMDsteps=$totMDsteps,tauT=$tauT,dt=$dt,thermostat=$thermostat,backupFreq=$backupFreq,heavyTrajFreq=$heavyTrajFreq,startfromzero=$startfromzero $scriptDIR/Thermalize.sbatch"
			sbatch --job-name=$nombre --export=exeDIR=$exeDIR,initConf=$initConf,Natoms=$Natoms,seed="${seed}",T=$T,totMDsteps=$totMDsteps,tauT=$tauT,dt=$dt,thermostat=$thermostat,backupFreq=$backupFreq,heavyTrajFreq=$heavyTrajFreq,startfromzero=$startfromzero $scriptDIR/Thermalize.sbatch
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
