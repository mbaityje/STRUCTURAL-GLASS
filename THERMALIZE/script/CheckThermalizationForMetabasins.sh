65;5401;1c#!/bin/bash
#
# Check thermalization of all the samples by calculating
# the self-intermediate scattering function, and tau.
#


#
# QUEUE PROPERTIES - SCHEDULER
#

if [ $(hostname) == 'banshee' ] || [ $(hostname) == 'PennPuter' ]; 
then readonly SCHEDULER='NONE'
elif [ -n $(hostname|grep talapas) ] || [ -n $(hostname|grep n[0-9][0-9][0-9]) ] || [ -n $(hostname|grep sn[0-9]) ]; 
then 
	readonly SCHEDULER='SLURM'
	queue=gpu
	simTime="1-00:00:00" #Days-HH:MM:SS
fi
readonly PROC_TAG="ctm"
readonly USERNAME=$(whoami)


#DIRECTORIES
rootDIR=$PWD/../..
echo rootDIR: $rootDIR
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
scriptDIR=$thermDIR/script
dataDIR=$thermDIR/data
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES
mkdir -p $workDIR
cd $workDIR
thermtimesFILE=$dataDIR/thermalizationtimes.txt

echo "Adesso mi trovo in $PWD"
readonly dt=0.0025
pot_type='KA'
pot_mode=${pot_mode:=shift}
thermostat=${thermostat:=NVE}

LISTAT=${1:-"5.0 2.0 1.0 0.8 0.7 0.6 0.49 0.46"}
LISTAN=${2:-"65"}
LISTASAM=${3:-"0"}

for T in $(echo $LISTAT)
do
	echo "T=$T"
	cd T$T
	
	for N in $(echo $LISTAN)
	do
		echo "N=$N"

		tautherm=$(awk -vT=$T -vN=$N '($1==T && $2==N){print $3}' $thermtimesFILE)
		nsteps=$(echo $tautherm/$dt|bc)

		cd N$N/$pot_mode

		if [ $LISTASAM == "all" ]; then
			LISTASAM=$(ls -d S[0-9]* | sed 's/S//')
		fi
		for ISAM in $(echo $LISTASAM)
		do
			cd S$ISAM
			echo "ISAM=$ISAM"

			thermalizeFile=thermalized.gsd
			if ! [ -f $thermalizeFile ]; then
				echo "$PWD/$thermalizeFile does not exist"
				cd ..
				continue
			fi

			if [ $SCHEDULER == "NONE" ]
			then
				thermostat=$thermostat pot_mode=$pot_mode pot_type=$pot_type rootDIR=$rootDIR bash $scriptDIR/SelfIntermediateScatteringFunction.sh $thermalizeFile 0 $nsteps $T $dt $tau_of_t
			elif [ $SCHEDULER == "SLURM" ]
			then
				nombre=N$N${PROC_TAG}T${T}i${ISAM}
				if [ 0 == `squeue -u$USERNAME -n $nombre|grep $USERNAME|wc -l` ]
				then
					echo "thermostat=$thermostat pot_mode=$pot_mode pot_type=$pot_type sbatch --job-name=$nombre -p$queue --export=filename=$thermalizeFile,iframe=0,nsteps=$nsteps,T=$T,dt=$dt,tau_of_t=$tau_of_t $scriptDIR/SelfIntermediateScatteringFunction.sh"
					sbatch --job-name=$nombre -p$queue --export=filename=$thermalizeFile,iframe=0,nsteps=$nsteps,T=$T,dt=$dt,tau_of_t=$tau_of_t,thermostat=$thermostat,pot_mode=$pot_mode,pot_type=$pot_type $scriptDIR/SelfIntermediateScatteringFunction.sh
				fi
			fi
			cd ..
		done #ISAM
		cd ../..
	done #N
	cd ..
done #T
