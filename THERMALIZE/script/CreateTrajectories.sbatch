#!/bin/bash
#SBATCH --ntasks=1

echo exeDIR: $exeDIR
echo filename: $filename
echo Natoms: $Natoms
echo seed: $seed
echo T: $T
echo pot_mode: $pot_mode
echo nsteps: $nsteps
echo tauT: $tauT
echo dt: $dt
echo thermostat: $thermostat
echo trajFreq: $trajFreq
echo backupFreq: $backupFreq
echo heavyTrajFreq: $heavyTrajFreq

python $exeDIR/ReadAndThermalize.py --user="$filename -N$Natoms -s$seed -T$T -t$nsteps --tau=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=0 --heavyTrajFreq=0 --trajFreq=$trajFreq --iframe=$iframe --addsteps -l$labelnew --startfromzero --pot_mode=$pot_mode --dumpacc"
echo "-------------------------------------"
echo Finished `basename "$0"`
echo "-------------------------------------"


