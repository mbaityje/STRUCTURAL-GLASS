#!/bin/bash
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
	simTime="0-02:00:00" #Days-HH:MM:SS
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
pot_mode='shift'
pot_type='KA'



for Tdir in T0.52 #T0.49 T0.48 T0.47 # T10.0 T0.52 T0.5 T5.0 T2.0 T1.0 T0.8 T0.7 T0.6 T0.55 #`ls -d T*|sort -r`
do
	T=`echo $Tdir | sed 's/^T//'`
	
	echo "-------------------------------------------"
	echo "|xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|"
	echo "-------------------------------------------"
	echo "T=$T"
	pwd
	cd $Tdir
	
	for Ndir in N65 #`ls -d N*`
	do
		N=`echo $Ndir | sed 's/^N//'`
		echo "N=$N"

		tautherm=$(awk -vT=$T -vN=$N '($1==T && $2==N){print $3}' $thermtimesFILE)
		nsteps=$(echo $tautherm/$dt|bc)

		cd $Ndir/$pot_mode
		for SAMdir in `ls -d S?`
		do
			ISAM=`echo $SAMdir | sed 's/^S//'`
			echo "ISAM=$ISAM"
			cd $SAMdir


			thermalizeFile=thermalized.gsd
			if ! [ -f $thermalizeFile ]; then
				echo "$PWD/$thermalizeFile does not exist"
				cd ..
				continue
			fi

			if [ $SCHEDULER == "NONE" ]
			then
				thermostat='NVT' pot_mode=$pot_mode pot_type=$pot_type rootDIR=$rootDIR bash $scriptDIR/SelfIntermediateScatteringFunction.sh $thermalizeFile 0 $nsteps $T $dt $tau_of_t
			elif [ $SCHEDULER == "SLURM" ]
			then
				nombre=N$N${PROC_TAG}T${T}i${ISAM}
				if [ 0 == `squeue -u$USERNAME -n $nombre|grep $USERNAME|wc -l` ]
				then
					thermostat='NVT' pot_mode=$pot_mode pot_type=$pot_type sbatch --job-name=$nombre -p$queue --export=filename=$thermalizeFile,iframe=0,nsteps=$nsteps,T=$T,dt=$dt,tau_of_t=$tau_of_t $scriptDIR/SelfIntermediateScatteringFunction.sh
				fi
			fi
			cd ..
		done #ISAM
		cd ../..
	done #N
	cd ..
done #T
