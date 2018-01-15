#!/bin/bash
#
# Take a trajectory starting from the frame $iframe.
# Minimize every single frame for $nframes.
# Save EIS, msdIS, nmovedIS.
#



T=${1:-0.43}
ISAM=${2:-0}
iframe=${3:-0}
nframes=${4:-10}
Natoms=65

readonly SYSTEM="PennPuter"
#readonly SYSTEM="Talapas"

#DIRECTORIES
rootDIR=$PWD/../..
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
scriptDIR=$thermDIR/script
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES

cd $workDIR/T$T/N$Natoms/S$ISAM
framesPerChunk=`python ~/STRUCTURAL-GLASS/UTILITIES/FindNFrames.py trajChunk0.gsd`
ichunk=`echo "$iframe/$framesPerChunk" |bc`

python $exeDIR/MinimizeSegment.py --user="trajChunk$ichunk.gsd --iframe=$iframe --nframes=$nframes"

