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


declare -A allM
allM=([5.0]=4 [2.0]=3  [1.0]=3  [0.8]=3  [0.7]=3  [0.6]=3  [0.55]=3  [0.52]=3  [0.5]=3 )

for T in 5.0 2.0 1.0 0.8 0.7 0.6 0.55
do
	M=${allM[$T]}
	for N in 1080
	do
		fr_diag=$(python $exeDIR/CalculateFriction.py $workDIR/T${T}/N${N}/Cd_NVE.txt --temperature=$T)
		fr_noise=$(python $exeDIR/CalculateFrictionJK.py $workDIR/T${T}/N${N}/noisecorrJK_NVT_M$M.txt $workDIR/T${T}/N${N}/noisecorrJK_NVT_M$M.npy --temperature=$T)
		echo "$T $N $fr_noise $fr_diag" >> $frictionFILE
	done
done

python ~/STRUCTURAL-GLASS/THERMALIZE/progs/CalculateFrictionJK.py noisecorrJK_NVT_M4.txt noisecorrJK_NVT_M4.npy -T5.0