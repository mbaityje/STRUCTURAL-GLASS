#!/bin/bash
#
# Take a trajectory starting from the frame $iframe.
# Minimize every single frame for $nframes.
# Save EIS, msdIS, nmovedIS.
#


readonly SYSTEM="PennPuter"
#DIRECTORIES
rootDIR=$PWD/../..
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
scriptDIR=$thermDIR/script
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES



TLIST=${1:-"0.8 0.7 0.55 0.52 0.49"}
SAMLIST=${2:-"0"}

Natoms=65

cd $workDIR
for T in $(echo $TLIST)
do
	for ISAM in $(echo $SAMLIST)
	do
		cd $workDIR/T$T/N$Natoms/shift/S$ISAM/chunksIS/
		iframe=0
		framesPerChunk=`python $utilDIR/FindNFrames.py trajChunk0.gsd`
		ichunk=`echo "$iframe/$framesPerChunk" |bc`	
		nframes=$framesPerChunk

		python $exeDIR/MinimizeSegment.py --user="trajChunk$ichunk.gsd --iframe=$iframe --nframes=$nframes" #--saveISgsd --saveThermalgsd

	done
done





exit




# T=${1:-0.6}
# ISAM=${2:-0}
# iframe=${3:-"default"}
# nframes=${4:-"default"}
# Natoms=65

# readonly SYSTEM="PennPuter"
# #readonly SYSTEM="talapas"

# #DIRECTORIES
# rootDIR=$PWD/../..
# thermDIR=$rootDIR/THERMALIZE
# exeDIR=$thermDIR/progs
# scriptDIR=$thermDIR/script
# workDIR=$rootDIR/OUTPUT
# utilDIR=$rootDIR/UTILITIES

# cd $workDIR/T$T/N$Natoms/shift/S$ISAM/chunksIS/
# framesPerChunk=`python $utilDIR/FindNFrames.py trajChunk0.gsd`
# ichunk=`echo "$iframe/$framesPerChunk" |bc`

# if [ $iframe  == "default" ]; then iframe=0; fi
# if [ $nframes == "default" ]; then nframes=$framesPerChunk; fi

# echo "nframes = $nframes"

# python $exeDIR/MinimizeSegment.py --user="trajChunk$ichunk.gsd --iframe=$iframe --nframes=$nframes" #--saveISgsd --saveThermalgsd

