#!/bin/bash
#SBATCH --job-name=therm
#SBATCH --ntasks=1
#SBATCH -p longgpu # partition (queue) 
#SBATCH --gres=gpu:1 
#SBATCH --time=0-02:00:00


echo exeDIR: $exeDIR
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

python $exeDIR/ReadAndThermalize.py --user="$filename -N$N -s$seed -T$T -t$nsteps --tau=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=0 --heavyTrajFreq=0 --trajFreq=$trajFreq --iframe=$iframe --addsteps -l$labelnew --startfromzero --pot_mode=$pot_mode --dumpacc"
echo "-------------------------------------"
echo Finished `basename "$0"`
echo "-------------------------------------"


