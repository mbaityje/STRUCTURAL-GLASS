#!/bin/bash
#
#A partire dalle configurazioni calde ne creo di piu` fredde

#
# QUEUE PROPERTIES - SCHEDULER
#

if [ $(hostname) == 'banshee' ] || [ $(hostname) == 'PennPuter' ]; 
then readonly SCHEDULER='NONE'
elif [ -n $(hostname|grep talapas) ] || [ -n $(hostname|grep n[0-9][0-9][0-9]) ] || [ -n $(hostname|grep sn[0-9]) ]; 
then 
	readonly SCHEDULER='SLURM'
	queue=gpu
	simTime="0-20:00:00" #Days-HH:MM:SS
fi
readonly PROC_TAG="rt"
readonly USERNAME=$(whoami)

#
# PHYSICAL AND NUMERICAL PARAMETERS
#
let nsamm1=$nsam-1
readonly hottestT=10.0
TLIST=${1:-"10.0 5.0 2.0 1.0 0.8 0.7 0.6 0.49 0.46"}
NLIST=${2:-"65"}
if [ -z $nsam       ]; then nsam=10; fi #default value of nsamples
if [ -z $pot_mode  ]; then pot_mode='shift'; fi #In these runs shift is the default mode
if [ -z $pot_type  ]; then pot_type='KA'; fi #In these runs Kob-Andersen is the default potential
if [ -z $dt         ]; then dt=0.0025 ; fi
if [ -z $backupFreq ]; then backupFreq=`echo 10/$dt|bc`; fi

echo backupFreq = $backupFreq

#
# DIRECTORIES
#
scriptDIR=$PWD
thermDIR=$PWD/..
rootDIR=$PWD/../..
exeDIR=$thermDIR/progs
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES
dataDIR=$thermDIR/data
mkdir -p $workDIR
cd $workDIR
thermtimesFILE=$dataDIR/thermalizationtimes.txt


#Range of temperatures
for T in $(echo $TLIST)
do
	mkdir -p T$T
	cd T$T
	echo "T = $T"

	for Natoms in $(echo $NLIST)
	do
		mkdir -p N$Natoms/$pot_mode
		cd N$Natoms/$pot_mode
		echo "Natoms = $Natoms, pot_mode: $pot_mode"
		
		inistateDIR=$workDIR/INITIAL-STATES/N$Natoms
		hottestTDIR=$workDIR/T$hottestT/N$Natoms/$pot_mode
		
		
		for isam in $(seq 0 $nsamm1)
		do
			echo "isam: $isam"
			mkdir -p S$isam
			cd S$isam

			seed=$(od -vAn -N4 -tu4 < /dev/urandom)
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
				totMDsteps=$(echo 14*$tauthermsteps|bc)
				thermostat=NVT
				heavyTrajFreq=$tauthermsteps
				tauT=0.1
				initConf=$hottestTDIR/S$isam/$thermConfName
			fi


			#If (partially) thermalized configuration already exists, then continue from there (to be able to make longer thermalization runs)
			echo Proposed initConf: $initConf
			if [ -e $thermConfName ] ;then
				initConf=$thermConfName
				donesteps=$(python $utilDIR/FindNsteps.py $initConf)
				if [ $donesteps -ge $totMDsteps ];
				then
					echo "This sample already has been thermalized, with $donesteps steps"
					cd ..; continue
				fi
				echo "Continuing run from $initConf"
				startfromzero=''
			else
				echo "Starting run from $initConf"
				startfromzero='--startfromzero'
			fi


			#It may happen that no initial configuration is available.
			#In that case, give a warning and skip sample
			if [ ! -f $initConf ]; 
			then
				echo "WARNING: $initconf does not exist"
				cd ..; continue
			fi

			#
			# Run Python script
			#
			if [ $SCHEDULER == "NONE" ]; 
			then
				echo `whoami`$USERNAME@`uname -n` > thermalized.time
				echo time "(python $exeDIR/ReadAndThermalize.py --user=\"$initConf -N$Natoms -s$seed -T$T -t$totMDsteps --tauT=$tauT --dt=$dt --pot_mode=$pot_mode --pot_type=$pot_type --thermostat=$thermostat --backupFreq=$backupFreq --heavyTrajFreq=$heavyTrajFreq $startfromzero\"  2>&1)" 2>>thermalized.time
				time (python $exeDIR/ReadAndThermalize.py --user="$initConf -N$Natoms -s$seed -T$T -t$totMDsteps --tauT=$tauT --dt=$dt --pot_mode=$pot_mode --pot_type=$pot_type --thermostat=$thermostat --backupFreq=$backupFreq --heavyTrajFreq=$heavyTrajFreq $startfromzero"  2>&1) 2>>thermalized.time
			elif [ $SCHEDULER == "SLURM" ]; 
			then
				nombre=N${Natoms}${PROC_TAG}T${T}i${isam}
				if [ 0 == `squeue -u$USERNAME -n $nombre|grep $USERNAME|wc -l` ]
				then
					sbatch --job-name=$nombre -p $queue --time=$simTime --export=exeDIR=$exeDIR,initConf=$initConf,Natoms=$Natoms,seed="${seed}",T=$T,pot_mode=$pot_mode,pot_type=$pot_type,totMDsteps=$totMDsteps,tauT=$tauT,dt=$dt,thermostat=$thermostat,backupFreq=$backupFreq,heavyTrajFreq=$heavyTrajFreq,startfromzero=$startfromzero,isam=$isam $scriptDIR/Thermalize.sbatch
				fi
			else
				echo "SCHEDULER = $SCHEDULER  not recognized"
				exit
			fi
			echo ""
			cd ..
		done #isam
		cd ../..
	done #Natoms
	cd ..
	wait
done #T
cd ..
