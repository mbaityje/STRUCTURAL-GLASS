#!/bin/bash
#SBATCH --job-name=therm
#SBATCH --ntasks=1
#SBATCH -p longgpu # partition (queue) 
#SBATCH --gres=gpu:1 
#SBATCH --time=0-06:00:00


echo exeDIR: $exeDIR
echo sample: $isam
echo initConf: $initConf
echo Natoms: $Natoms
echo seed: $seed
echo T: $T
echo potential: $potential
echo totMDsteps: $totMDsteps
echo tauT: $tauT
echo dt: $dt
echo thermostat: $thermostat
echo backupFreq: $backupFreq
echo heavyTrajFreq: $heavyTrajFreq
echo startfromzero: $startfromzero


python $exeDIR/ReadAndThermalize.py --user="$initConf -N$Natoms -s$seed -T$T -t$totMDsteps --tauT=$tauT --pot_mode=$potential --dt=$dt --thermostat=$thermostat --backupFreq=$backupFreq --heavyTrajFreq=$heavyTrajFreq startfromzero" 

echo "-------------------------------------"
echo Finished `basename "$0"`
echo "-------------------------------------"


