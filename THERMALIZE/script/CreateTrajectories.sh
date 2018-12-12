#!/bin/bash
#
# Takes as input some configurations which are assumed to be well-thermalized.
# For each of them, runs and saves a series of trajectories.
# 
#

module switch anaconda3 anaconda3/4.4.0 #The version of cudatoolkit in future builds is incompatible with hoomd

#Relevant directories
#
readonly PROC_TAG=crtr
readonly USERNAME=$(whoami)
echo We are at: $(pwd)


if [ `hostname` == "PennPuter" ] || [ `hostname` == "banshee" ];
then SYSTEM="PennPuter";
elif [ `hostname` == "tango" ];
then SYSTEM="PennPuter";
else SYSTEM="talapas"; fi
echo SYSTEM = $SYSTEM

#DIRECTORIES
if [ $SYSTEM == "talapas" ];
then 
    rootDIR=/home/mbaity/STRUCTURAL-GLASS/
    if [ -z $simTime ]; then simTime="0-06:00:00"; fi
    if [ -z $queue   ]; then   queue=gpu; fi
    if [[ $queue == *"gpu"* ]];
    then
	gres="--gres=gpu:1";
	echo "Reserving a GPU";
    else
	echo "Reserving only CPU (no GPU)"
    fi
else 
    rootDIR=$PWD/../..
    echo "Setting rootDIR: $rootDIR"
fi
thermDIR=$rootDIR/THERMALIZE
scriptDIR=$thermDIR/script
exeDIR=$thermDIR/progs
dataDIR=$thermDIR/data
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES
cd $workDIR
thermtimesFILE=$dataDIR/thermalizationtimes.txt

#Parameters are taken from command line, though I set some default value, mainly for code developing and debugging
ntraj=${1:-2} #Number of trajectories we want, for each choice of the parameters
let ntrajm1=$ntraj-1
LISTAT=${2:-"5.0"}
LISTAN=${3:-"1080"}
LISTAISAM=${4:-"0"}

echo "ntraj     = $ntraj"
echo "lista T   : $LISTAT"
echo "lista N   : $LISTAN"
echo "lista isam: $LISTAISAM"

if [ -z $thermostat ]; then thermostat='NVT'; fi
if [ -z $pot_mode ]; then pot_mode='xplor'; fi


readonly tauT=0.1
readonly dt=0.0025
trajFreq=-1000 #A negative trajFreq means we sample |trajFreq| quasi-logarithmically distributed bins

echo "thermostat: $thermostat"
echo "potential mode: $pot_mode"
echo "dt: $dt"


for T in $(echo $LISTAT)
do
	for N in $(echo $LISTAN)
	do
		#My simulations of bigger systems always have an xplor potential
		if [ $N -gt 500 ]; then pot_mode='xplor'; echo 'Script decides that potential is xplor'; fi

		tautherm=$(awk -vT=$T -vN=$N '($1==T && $2==N){print $3}' $thermtimesFILE)
		nsteps=$(echo 2*$tautherm/$dt|bc) # I put an extra factor 1.5 inside it so that there is a small bit with no signal

		for isam in $(echo $LISTAISAM)
		do
			directorio=$workDIR/T$T/N$N/S$isam/trajectories/
			mkdir -p $directorio
			cd $directorio
			ln -s ../thermalized.gsd thermalized.gsd 2> /dev/null #Silence stderr otherwise it bitches everytime the link alread exists

			nextitraj=`ls ${thermostat}${potential}*.gsd 2> /dev/null| wc -l`
			echo nextitraj: $nextitraj
			for itraj in $(seq $nextitraj $ntrajm1)
			do
				echo ""; echo T=$T N=$N isam=$isam itraj=$itraj
				labelold=${thermostat}${potential}$(echo $itraj-1|bc) #label for the read configuration
				labelnew=${thermostat}${potential}${itraj} #label for the configuration that will be generated
				if [ ${itraj} == 0 ]; then filename=thermalized.gsd; else filename=${labelold}.gsd; fi
				echo filename: $filename
				let iframe=$(python $utilDIR/FindNFrames.py $filename)-1

				if [ $SYSTEM == "PennPuter" ]; 
				then
					python $exeDIR/ReadAndThermalize.py --user="$filename -N$N -s0 -T$T -t$nsteps --tau=$tauT --pot_mode=$pot_mode --dt=$dt --thermostat=$thermostat --backupFreq=10000 --heavyTrajFreq=0 --trajFreq=$trajFreq --iframe=$iframe --addsteps -l$labelnew --startfromzero --dumpacc"


				elif [ $SYSTEM == "talapas" ]; 
				then
					nombre=N${N}${PROC_TAG}T${T}i${isam}
					if [ 0 == `squeue -u$USERNAME -n $nombre|grep $USERNAME|wc -l` ]
					then
						sbatch $gres --time=$simTime -p$queue --job-name=$nombre --export=all,exeDIR=$exeDIR,filename=$filename,Natoms=$N,seed=0,T=$T,nsteps=$nsteps,tauT=$tauT,pot_mode=$pot_mode,dt=$dt,thermostat=$thermostat,backupFreq=10000,heavyTrajFreq=0,trajFreq=$trajFreq,iframe=$iframe,addsteps=$addsteps,labelnew=$labelnew,dumpacc=$dumpacc $scriptDIR/CreateTrajectories.sbatch
					fi
				else
					echo "SYSTEM=$SYSTEM not recognized"
					exit
				fi



			done
		done
	done
done
