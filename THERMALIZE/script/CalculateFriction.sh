#!/bin/bash


SYSTEM="PennPuter"

#DIRECTORIES
if [ $SYSTEM == "talapas" ];
then rootDIR=/home/mbaity/STRUCTURAL-GLASS/
else 
	rootDIR=$PWD/../..
	echo "Setting rootDIR: $rootDIR"
fi
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
dataDIR=$thermDIR/data
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES

frictionFILE=$dataDIR/friction.txt
echo "T N friction_noise friction_diag" > $frictionFILE

for T in 5.0 2.0 1.0 0.6
do
	for N in 1080
	do
		fr_diag=$(python $exeDIR/CalculateFriction.py $workDIR/T${T}/N${N}/Cd_NVE.txt --temperature=$T)
		fr_noise=$(python $exeDIR/CalculateFriction.py $workDIR/T${T}/N${N}/noisecorr_NVE.txt --temperature=$T)
		echo "$T $N $fr_noise $fr_diag" >> $frictionFILE
	done
done

