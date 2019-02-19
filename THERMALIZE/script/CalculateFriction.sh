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

LISTAT=${1:-"5.0"}
LISTAN=${2:-"1080"}
LISTATHERMOSTAT=${3:-"NVT"}
if [ $showplots ]; then showplots="--showplots"; fi

frictionFILE=$dataDIR/friction.txt
echo "T N friction_noise friction_diag" > $frictionFILE


declare -A allM
#allM=([5.0]=4 [2.0]=3  [1.0]=3  [0.8]=3  [0.7]=3  [0.6]=3  [0.55]=3  [0.52]=3  [0.5]=3 ) #For JK

for T in $(echo $LISTAT)
do
	echo "T = $T"
	#	M=${allM[$T]} #for JK
	M=3
	for N in $(echo $LISTAN)
	do
		for thermostat in $(echo $LISTATHERMOSTAT)
		do	
			echo "Only mean values"
			python $exeDIR/CalculateFrictionNoise.py $workDIR/T${T}/N${N}/noisecorr_NVT_combine_M$M.txt --temperature=$T --thermostat=$thermostat
			python $exeDIR/CalculateFrictionDiag.py  $workDIR/T${T}/N${N}/Cd_NVT.txt --temperature=$T --thermostat=$thermostat

			echo "Proper jackknife"
			python $exeDIR/CalculateFrictionNoiseJK.py $workDIR/T${T}/N${N}/noisecorrJK_NVT_M$M.txt $workDIR/T${T}/N${N}/noisecorrJK_NVT_M$M.npy --temperature=$T --thermostat=$thermostat $showplots

		# fr_diag=$(python $exeDIR/CalculateFriction.py $workDIR/T${T}/N${N}/Cd_NVT.txt --temperature=$T)
		# fr_noise=$(python $exeDIR/CalculateFriction.py $workDIR/T${T}/N${N}/noisecorr_NVT_combine_M$M.txt --temperature=$T)
		# #fr_noise=$(python $exeDIR/CalculateFrictionJK.py $workDIR/T${T}/N${N}/noisecorrJK_NVT_M$M.txt $workDIR/T${T}/N${N}/noisecorrJK_NVT_M$M.npy --temperature=$T)
		# echo "$T $N $fr_noise $fr_diag" >> $frictionFILE
		done
	done
done

