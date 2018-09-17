#!/bin/bash
#SBATCH --ntasks=1
#SBATCH -p longgpu # partition (queue)
#SBATCH --gres=gpu:1
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


if [ `hostname` == "PennPuter" ];
then SYSTEM="PennPuter";
elif [ `hostname` == "tango" ];
then SYSTEM="PennPuter";
else SYSTEM="talapas"; fi

#DIRECTORIES
if [ $SYSTEM == "talapas" ];
then rootDIR=/home/mbaity/STRUCTURAL-GLASS/
else 
	rootDIR=$PWD/../../../..
	echo "Setting rootDIR: $rootDIR"
fi
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES

#
# Command line input
#
# If the script was launched with sbatch, the following variables are automatically
# assigned. If it was launched with bash, they are assigned via $1, $2, ...
filename=${filename:-$1}
iframe=${iframe:-$2}
nsteps=${nsteps:-$3}
T=${T:-$4}
dt=${dt:-$5}
tau_of_t=${tau_of_t:-$6}

#Set default values for the last two arguments [i.e. if they were unassigned give them default value]
dt=${dt:-0.0025}
tau_of_t=${tau_of_t:-0} #flag that tells us if we are seeking to calculate only one tau (tau_of_t=1) or if we want to calculate tau two times, spaced by a gap (tau_of_t=0)

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
maxFrames=1000 #The (first) trajectory we construct has at most 1000 frames
ratio=`echo "$nsteps/$maxFrames" | bc`
trajFreq=-1000 #A negative trajFreq means we sample -trajFreq logarithmically distributed bins

#
#Some checks to make sure that the input is good
#
#Number of arguments
if [ $SYSTEM == "PennPuter" ] && [ $# -ne 4 ] && [ $# -ne 5 ] && [ $# -ne 6 ]; then
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
pot_mode='xplor'
label="_ifr${iframe}_${pot_mode}"
rm -f trajectory${label}.gsd
echo python $exeDIR/ReadAndThermalize.py --user=\"$filename -N$Natoms -s0 -T$T -t$nsteps --tau=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=0 --heavyTrajFreq=0 --iframe=$iframe --trajFreq=$trajFreq --addsteps --startfromzero --pot_mode=$pot_mode\"
python $exeDIR/ReadAndThermalize.py --user="$filename -N$Natoms -s0 -T$T -t$nsteps --tau=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=0 --heavyTrajFreq=0 --iframe=$iframe --trajFreq=$trajFreq --addsteps -l$label --startfromzero --pot_mode=$pot_mode"

#
# Calculate Fk(t), MSD and tau for the first time
#
echo "|--> Calculating first Fk(t)..."
echo "python $exeDIR/SelfIntermediateScatteringFunction.py  trajectory${label}.gsd --dt=$dt --every_forMemory=1 -l${label}"
python $exeDIR/SelfIntermediateScatteringFunction.py  trajectory${label}.gsd --dt=$dt --every_forMemory=1 -l${label}


if [ 0 -eq $tau_of_t ]
then
	#
	# Now run a gap of 20tau without saving any backup
	#
	echo "\n|--> Running gap of 20 tau..."
	filenamegap=$label.gsd #We read from the output of the previous simulation, which is $label.gsd
	nstepsgap=`echo 20*${nsteps} | bc`
	labelgap="_gap_${pot_mode}"
	python $exeDIR/ReadAndThermalize.py --user="$filenamegap -N$Natoms -s0 -T$T -t$nstepsgap --tau=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=0 --heavyTrajFreq=0 --trajFreq=0 --iframe=$iframe --addsteps -l$labelgap --startfromzero --pot_mode=$pot_mode"

	#
	# Run 3tau saving the trajectory
	#
	echo "\n|--> Running new trajectory for a new measurement of Fk(t)..."
	filenameaftergap=$labelgap.gsd #We read from the output of the previous simulation, which is $labelgap.gsd
	labelaftergap="_aftergap_${pot_mode}"
	nstepsaftergap=$nsteps
	
	rm -f trajectory${labelaftergap}.gsd
	python $exeDIR/ReadAndThermalize.py --user="$filenameaftergap -N$Natoms -s0 -T$T -t$nstepsaftergap --tau=$tauT --dt=$dt --thermostat=$thermostat --backupFreq=0 --heavyTrajFreq=0 --iframe=$iframe --trajFreq=$trajFreq --addsteps -l$labelaftergap --pot_mode=$pot_mode"
	
	
	#
	# Calculate Fk(t), MSD and tau for the second time
	#
	echo "\n|--> Calculating Fk(t) again..."
	python $exeDIR/SelfIntermediateScatteringFunction.py  trajectory${labelaftergap}.gsd --dt=$dt --every_forMemory=1 -l${labelaftergap}

else
	echo "tau_of_t=1 not implemented"
	exit
fi

echo "+++++++++++++++++++++++++++++++++++++++"
