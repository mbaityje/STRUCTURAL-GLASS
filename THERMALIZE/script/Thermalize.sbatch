#!/bin/bash
#SBATCH --job-name=therm
#SBATCH --ntasks=1
#SBATCH -p longgpu # partition (queue) 
#SBATCH --gres=gpu:1 


echo exeDIR: $exeDIR
echo initConf: $initConf
echo Natoms: $Natoms
echo seed: $seed
echo T: $T
echo totMDsteps: $totMDsteps
echo tauT: $tauT
echo dt: $dt
echo thermostat: $thermostat
echo backupFreq: $backupFreq
echo heavyTrajFreq: $heavyTrajFreq

module list
ls /packages/anaconda2/4.4.0/bin/python
which python
python --version


python $exeDIR/ReadAndThermalize.py --user="$initConf -N$Natoms -s$seed -T$T -t$totMDsteps --tauT=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=$backupFreq --heavyTrajFreq=$heavyTrajFreq" 

echo "-------------------------------------"
echo Finished `basename "$0"`
echo "-------------------------------------"


