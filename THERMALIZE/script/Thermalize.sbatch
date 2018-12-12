#!/bin/bash
#SBATCH --ntasks=1

echo exeDIR: $exeDIR
echo sample: $isam
echo initConf: $initConf
echo Natoms: $Natoms
echo seed: $seed
echo T: $T
echo pot_mode: $pot_mode
echo totMDsteps: $totMDsteps
echo tauT: $tauT
echo dt: $dt
echo thermostat: $thermostat
echo backupFreq: $backupFreq
echo heavyTrajFreq: $heavyTrajFreq
echo startfromzero: $startfromzero


python $exeDIR/ReadAndThermalize.py --user="$initConf -N$Natoms -s$seed -T$T -t$totMDsteps --tauT=$tauT --pot_mode=$pot_mode --dt=$dt --thermostat=$thermostat --backupFreq=$backupFreq --heavyTrajFreq=$heavyTrajFreq $startfromzero" 

echo "-------------------------------------"
echo Finished `basename "$0"`
echo "-------------------------------------"


