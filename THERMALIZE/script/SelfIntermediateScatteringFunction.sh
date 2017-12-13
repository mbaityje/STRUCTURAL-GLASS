#!/bin/bash
#
# Given a gsd file and a target frame:
# - Runs a trajectory starting from that configuration
# - Calculates the self-intermediate scattering function
# - Calculates the relaxation time from the self-intermediate scattering function
# If tau_of_t==0
# -runs again for 20tau
# -runs again a trajectory
# -calculates again the observables
# -compares the result between the runs to see whether we thermalized
# 

#Relevant directories
#
echo "We are at"
pwd


#SYSTEM="PennPuter"
SYSTEM="Talapas"
if [ $SYSTEM == "Talapas" ];
then rootDIR=/home/mbaity/STRUCTURAL-GLASS/
else 
    echo "Implement rootDIR for SYSTEM!=talapas"
    exit
fi
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES

#
#Command line input
#
filename=$1
iframe=$2
nsteps=$3
T=$4
dt=${5:-0.0025} #dt is taken from command line. Otherwise it is 0.0025
tau_of_t=${6:-0} #flag that tells us if we are seeking to calculate only one tau (tau_of_t=1) or if we want to calculate tau two times, spaced by a gap (tau_of_t=0)

echo "SelfIntermediateScatteringFunction.sh:"
echo "filename: $filename"
echo "iframe: $iframe"
echo "nsteps: $nsteps"
echo "T = $T"
echo "dt = $dt"
echo "tau_of_t: $tau_of_t"

#Some hardcoded parameters that I might decide to put as command-line input
readonly thermostat='NVE'
readonly tauT=0.1
readonly Natoms=65
maxFrames=1000 #The first trajectory we construct has at most 1000 frames
ratio=`echo "$nsteps/$maxFrames" | bc`
trajFreq=$((ratio<1?1:ratio))
#trajFreq=1


#
#Some checks to make sure that the input is good
#
#Number of arguments
if [ $# -ne 4 ] && [ $# -ne 5 ] && [ $# -ne 6 ]; then
	echo "Wrong number of parameters ($#)"
	echo "Launch as: $0 <filename> <iframe> <nsteps> <T>"
	exit
fi
#Filename
if ! [ -e $filename ]; then	
	echo " \"$filename\" is not a readable file"; 
	echo "Launch as: $0 <filename> <iframe> <nsteps> <T>"
	exit
fi
#Frame index must be a positive integer
if ! [[ $iframe =~ ^[0-9]+$ ]]; then
	echo "iframe (=$iframe) must be a non-negative integer"
	echo "Launch as: $0 <filename> <iframe> <nsteps> <T>"
	exit
else
	#Frame index must be less than the number of frames
	Nframes=`python $utilDIR/FindNFrames.py $filename`
	if [ $iframe -ge $Nframes ]; then
		echo "iframe ($iframe) cannot be larger than the total number of frames in the file ($Nframes)"
		echo "Launch as: $0 <filename> <iframe> <nsteps> <T>"
	fi
fi
#nsteps must be a positive integer
if ! [[ $iframe =~ ^[0-9]+$ ]]; then
	echo "nsteps (=$nsteps) must be a non-negative integer"
	echo "Launch as: $0 <filename> <iframe> <nsteps> <T>"
	exit
fi
#T must be a non negative number
if ! [[ $T =~ ^[0-9]*\.?[0-9]*$ ]]; then
	echo "T (=$T) must be a non-negative integer"
	echo "Launch as: $0 <filename> <iframe> <nsteps> <T>"
	exit
fi


#
#Now we run a simulation of the chosen length, in which the configuration is dumped often
#
echo "|--> Creating first trajectory..."
addsteps='True'
label="_ifr$iframe"
rm -f trajectory${label}.gsd
echo python $exeDIR/ReadAndThermalize.py --user=\"$filename -N$Natoms -s0 -T$T -t$nsteps --tau=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=0 --heavyTrajFreq=0 --iframe=$iframe --trajFreq=$trajFreq --addsteps=$addsteps\"
python $exeDIR/ReadAndThermalize.py --user="$filename -N$Natoms -s0 -T$T -t$nsteps --tau=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=0 --heavyTrajFreq=0 --iframe=$iframe --trajFreq=$trajFreq --addsteps=$addsteps -l$label"

#
# Calculate Fk(t), MSD and tau for the first time
#
echo "|--> Calculating first Fk(t)..."
echo "python $exeDIR/SelfIntermediateScatteringFunction.py  trajectory${label}.gsd --dt=$dt --every_forMemory=1 -l${label}"
python $exeDIR/SelfIntermediateScatteringFunction.py  trajectory${label}.gsd --dt=$dt --every_forMemory=1 -l${label}

tauFkt_file=tau$label.txt
if [ -f $tauFkt_file ];
then 
    tauFkt=`tail -1 $tauFkt_file`
    echo "TAU is $tauFkt"
else
    echo "We were not able to calculate TAU. We will go on, and calculate all the other quantities. If the problem is only the fluctuations, you will be able to calculate TAU once you averaged over the samples."
    echo "We assume TAU=$nsteps"
    tauFkt=$nsteps
fi



if [ 0 -eq $tau_of_t ]
then
    #
    # Now run a gap of 20tau without saving any backup
    #
    echo "|--> Running gap of 20 tau..."
    addsteps='True'
    filenamegap=$label.gsd #We read from the output of the previous simulation, which is $label.gsd
    nstepsgap=`echo 20*${tauFkt}/$dt | bc`
    labelgap='_gap'
    python $exeDIR/ReadAndThermalize.py --user="$filenamegap -N$Natoms -s0 -T$T -t$nstepsgap --tau=$tauT --dt=$dt --thermostat=$thermostat \
											--backupFreq=0 --heavyTrajFreq=0 --trajFreq=0 --iframe=$iframe --addsteps=$addsteps -l$labelgap"

    #
    # Run 3tau saving the trajectory
    #
    echo "|--> Running trajectory of 3 tau..."
    addsteps='True'
    filenameaftergap=$labelgap.gsd #We read from the output of the previous simulation, which is $labelgap.gsd
    labelaftergap='_aftergap'
    nstepsaftergap=`echo 3*$tauFkt/$dt | bc`
    
    rm -f trajectory${labelaftergap}.gsd
    python $exeDIR/ReadAndThermalize.py --user="$filenameaftergap -N$Natoms -s0 -T$T -t$nstepsaftergap --tau=$tauT --dt=$dt --thermostat=$thermostat \
												--backupFreq=0 --heavyTrajFreq=0 --iframe=$iframe --trajFreq=$trajFreq --addsteps=$addsteps -l$labelaftergap"
    
    
    #
    # Calculate Fk(t), MSD and tau for the second time
    #
    echo "|--> Calculating Fk(t) again..."
    python $exeDIR/SelfIntermediateScatteringFunction.py  trajectory${labelaftergap}.gsd --dt=$dt --every_forMemory=1 -l${labelaftergap}
    tauFkt2_file=tau$labelaftergap.txt
    if [ -f $tauFkt2_file ];
    then 
	tauFkt2=`tail -1 $tauFkt2_file`
    fi
    echo "TAU1 is $tauFkt"
    echo "TAU2 is $tauFkt2"
    
    #
    # State if the sample is thermalized
    #
    #I say it is thermalized if they are within 5%
    rtol=0.10
    rdiff=`echo "($tauFkt-$tauFkt2)/$tauFkt2" | bc -lq`
    absrdiff=`echo "sqrt($rdiff*$rdiff)" | bc -lq`
    
    #Create a dummy file to quickly check if the directory is thermalized
    rm -f YES_thermalized NOT_thermalized
    if [ 1 == $( echo "$absrdiff <= $rtol"|bc ) ]; 
    then
	echo "Thermalized"
	date > YES_thermalized
	echo "|tau1-tau2|/tau2 = $absrdiff <= $rtol" > YES_thermalized
	
    else
	echo "Not thermalized"
	date > NOT_thermalized
	echo "|tau1-tau2|/tau2 = $absrdiff > $rtol" > NOT_thermalized
    fi
else
    echo "tau_of_t not really implemented"
    exit
fi

echo "+++++++++++++++++++++++++++++++++++++++"
