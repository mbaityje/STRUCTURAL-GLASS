#!/bin/bash
#Process tag is used for the queue scheduler
readonly PROC_TAG="im"

if [ `hostname` == "PennPuter" ] || [ `hostname` == "tango" ] || [ `hostname` == "banshee" ];
then
	SCHEDULER="NONE";
else
	SCHEDULER="SLURM";
	queue=gpu
	simTime="1-00:00:00" #Days-HH:MM:SS
fi

#DIRECTORIES
rootDIR=$PWD/../..
thermDIR=$rootDIR/THERMALIZE
scriptDIR=$thermDIR/script
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES
exeDIR=$thermDIR/progs


#Parameters
TLIST=${1:-"0.6"}
pot_mode=${pot_mode:-"shift"}
SAMLIST=${SAMLIST:-"0"}
thres=${thres:-1e-5}
maxbasins=${maxbasins:-3000}
readonly Natoms=65
readonly dt=0.0025

if [ $showplots ]; then showplots="--showplots"; else showplots=""; fi

tauNAME=tauMB.txt
qNAME=qMB.txt
eNAME=EMB.txt
eridgeNAME=EridgeMB.txt


cd $workDIR
for T in $(echo $TLIST)
do
	cd T$T/N$Natoms/$pot_mode
	rm  $tauNAME $qNAME $eNAME $eridgeNAME
	for ISAM in $(echo $SAMLIST) #$(ls -d S* | sed 's/S//') #ISAM=any sample that has been simulated
	do

		echo -e T = $T      ISAM = $ISAM
		cd S$ISAM/chunksIS
		L=$(python $utilDIR/FindL.py ../thermalized.gsd) #3.78364777565

		python $exeDIR/IdentifyMetabasins.py ./ -L$L --thres=$thres --maxbasins=$maxbasins $showplots

		# Create a single file with all the samples inside it
		cat tauMB.txt >> ../../$tauNAME
		cat qMB.txt >> ../../$qNAME
		cat EMB.txt >> ../../$eNAME
		cat EridgeMB.txt >> ../../$eridgeNAME

		cd ../..
	done
	pwd
	# alltau=$(awk '{print $0}' S*/chunksIS/tauMB.txt)
	# echo $alltau > ~/Borrar/tau.txt
	cd ../../..
done
pwd
echo $workDIR

