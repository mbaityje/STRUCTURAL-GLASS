#!/bin/bash
#
# Check thermalization of all the samples by calculating
# the self-intermediate scattering function, and tau.
#

readonly PROC_TAG="ct"
readonly USERNAME=`whoami`
queue=gpu
#queue=short


#Script Options
readonly tau_of_t=0 #1: calculate Fkt on all the heavyTraj, 0: calculate Fkt on only the last configuration

if [ `hostname` == "PennPuter" ] || [ `hostname` == "banshee" ];
then SYSTEM="PennPuter";
elif [ `hostname` == "tango" ];
then SYSTEM="PennPuter";
else SYSTEM="talapas"; fi

#DIRECTORIES
rootDIR=$PWD/../..
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



for Tdir in T0.52 #T0.7 #T0.55 #T0.8 T5.0 T2.0 T1.0 T0.6 #`ls -d T*|sort -r`
do
	T=`echo $Tdir | sed 's/^T//'`
	
	echo "-------------------------------------------"
	echo "|xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|"
	echo "-------------------------------------------"
	echo "T=$T"
	pwd
	cd $Tdir
	
	for Ndir in N1080 #`ls -d N*`
	do
		N=`echo $Ndir | sed 's/^N//'`
		echo "N=$N"

		tautherm=$(awk -vT=$T -vN=$N '($1==T && $2==N){print $3}' $thermtimesFILE)
		nsteps=$(echo $tautherm/$dt|bc)

		cd $Ndir
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

			if [ $SYSTEM == "PennPuter" ]
			then
				pot_mode='xplor' bash $scriptDIR/SelfIntermediateScatteringFunction.sh $thermalizeFile 0 $nsteps $T $dt $tau_of_t
			elif [ $SYSTEM == "talapas" ]
			then
				nombre=N$N${PROC_TAG}T${T}i${ISAM}
				if [ 0 == `squeue -u$USERNAME -n $nombre|grep $USERNAME|wc -l` ]
				then
					sbatch --job-name=$nombre -p$queue --export=filename=$thermalizeFile,iframe=0,nsteps=$nsteps,T=$T,dt=$dt,tau_of_t=$tau_of_t $scriptDIR/SelfIntermediateScatteringFunction.sh
				fi
			fi
			cd ..
		done #ISAM
		cd ..
	done #N
	cd ..
done #T
