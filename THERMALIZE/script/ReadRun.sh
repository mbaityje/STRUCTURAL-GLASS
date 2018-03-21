### #!/bin/bash

#This script launches ReadRun.py for the given set of samples


for isam in $(echo $batch)
do
	echo "sam: $isam"
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
	echo `whoami`$USERNAME@`uname -n` > thermalized.time
	echo python $exeDIR/ReadRun.py --user=\"$initConf -p$paramsFILE\"
	time (python $exeDIR/ReadRun.py --user="$initConf -p$paramsFILE" 2>&1) 2>>thermalized.time

	cd $workDIR
done

