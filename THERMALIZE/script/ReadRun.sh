#!/bin/bash
#SBATCH --job-name=readrun
#SBATCH --ntasks=1
#SBATCH -p gpu # partition (queue)
#SBATCH --gres=gpu:1
#SBATCH --time=0-01:00:00
##SBATCH -o ~/STRUCTURAL-GLASS/LOGS/job$(SLURM_JOB_ID).txt

################################################################
#                                                              #
# This script launches ReadRun.py for the given set of samples #
# Environment variables are passed through RunAll.sh           #
#                                                              #
################################################################


echo "-- Comincia readRun.sh --"
source /etc/profile.d/modules.sh
module load anaconda3/4.4.0

for isam in $(echo $batch)
do
	echo "sam: $isam"
	samOUT=ReadRun$$.log
	mkdir -p $workDIR/T$T/N$Natoms/S$isam
	cd $workDIR/T$T/N$Natoms/S$isam

	#If a partially thermalized state exists, start from there.
	#Otherwise, start from T=10 (unless this run is at T=10, then start from crystal)
	thermConfName=thermalized.gsd
	if [ -e $thermConfName ] ;then initConf=$thermConfName
	else
		if [ $T == "$hottestT" ]; then initConf=$inistateDIR/initIS.gsd
		else initConf=$hottestTDIR/S$isam/$thermConfName; fi
		if [ ! -f $initConf ]; then
			echo "WARNING: $initconf does not exist"
			cd $workDIR; continue
		fi
	fi

	#Execution
	pwd
	echo "File di output: $PWD/$samOUT"
	echo "SYSTEM: $SYSTEM"
	if [ $SYSTEM == 'PennPuter' ]; then
	    echo "SYSTEM=PennPuter. Launching minibatch serially."
	    echo "python $exeDIR/ReadRun.py --user=\"$initConf -p$paramsFILE\""
	    python $exeDIR/ReadRun.py --user="$initConf -p$paramsFILE"
	elif [ $SYSTEM == 'talapas' ];  then
		echo "python $exeDIR/ReadRun.py --user=\"$initConf -p$paramsFILE\"" >> $samOUT
		python $exeDIR/ReadRun.py --user="$initConf -p$paramsFILE" >> $samOUT &
		echo "python script finished" >> $samOUT
	else 
		echo "SYSTEM=$SYSTEM is not a recognised option. Choose between PennPuter and talapas"
		exit
	fi
	cd $workDIR
done

echo "Waiting..."
wait
