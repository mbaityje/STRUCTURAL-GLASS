#!/bin/bash

#
#Command line input. Default values are for testing purposes. 
#
T=${1}
dt=${2}
deltaE=${3}
ttot=${4}
doTS=${5:-0}

#
#Some checks to make sure that the input is good
#
#Number of arguments
if [ $# -ne 4 ] && [ $# -ne 5 ]; then
	echo "Wrong number of parameters ($#)"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doTS>]"
	echo "T: temperature"
	echo "dt: MD integration step"
	echo "deltaE: energy difference to state that two neighboring IS are different"
	echo "ttot: total length of the trajectory"
	echo "doTS: 0 (default): Calculate the trajectory of the IS. 1: Calculate both IS and TS."
	exit
fi
#T must be a non negative number
if ! [[ $T =~ ^[0-9]*\.?[0-9]*$ ]]; then
	echo "T (=$T) must be a non-negative number"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doTS>]"
	exit
fi
#dt must be 0<dt<0.01
if ! [[ $dt =~ ^[0-9]*\.?[0-9]*$ ]]; then
	echo "dt (=$deltaE) must be a non-negative number"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doTS>]"
	exit
elif [ 1 -eq $(echo "$dt > 0.01" | bc -lq)  ]; then
	echo "dt = $dt >= 0.01, is very big, you will likely have particles falling out of the box"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doTS>]"
	exit
fi
#deltaE must be a non negative number
if ! [[ $deltaE =~ ^[0-9]*\.?[0-9]*$ ]]; then
	echo "deltaE (=$deltaE) must be a non-negative integer"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doTS>]"
	exit
fi
#ttot must be a positive integer
if ! [[ $ttot =~ ^[0-9]+$ ]]; then
	echo "ttot (=$ttot) must be a positive integer"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doTS>]"
	exit
fi
#doTS must be 0 or 1
if ! [ $doTS -eq 0 -o $doTS -eq 1 ]; then
	echo "doTS ($doTS) must be 0 or 1."
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doTS>]"
	exit
fi


#Hard-coded parameters
readonly tauT=0.1
readonly tchunk=`echo 10^5|bc`
nchunks=`echo $ttot/$tchunk|bc`
let nchunksm1=$nchunks-1

#DIRECTORIES
rootDIR=$PWD/../../../..
thermDIR=$rootDIR/THERMALIZE
exeDIR=$thermDIR/progs
scriptDIR=$thermDIR/script
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES


#--------------------------------------------#
# Here starts the actual bisection by chunks #
#--------------------------------------------#

# Filenames
filename=thermalized.gsd
elistFILE=elist.txt

#Initial time and chunk (if a previous run was interrupted we restart from where it finished)
tlast=0
firstChunk=0
if [ -f $elistFILE ];
then
	if [ `wc -l $elistFILE |cut -f1 -d" "` -gt 0 ]
	then
	    tlast=`tail -1 $elistFILE|cut -f1 -d" "`;
	    firstChunk=`echo "($tlast+1)/$tchunk"|bc`
	fi
fi


#Iterate over the chunks
for ichunk in $(seq $firstChunk $nchunksm1)
do
	echo "+++ ichunk = $ichunk +++"

	#First generate the thermal trajectory of the cunk
	python $exeDIR/CreateChunk.py --user="$filename --ichunk=$ichunk --tchunk=$tchunk --dt=$dt --temperature=$T --tauT=$tauT"

	#Now that we have the new chunk, we can delete the previous one
	#Uncomment the next line if we don't want to keep the thermal trajectory
	#if [ $ichunk -ge 1 ]; then rm -f restartChunk`expr $ichunk - 1`.gsd; fi
	
	#Perform IS bisection on the chunk
	echo "python $exeDIR/BisectChunk.py --user=\"trajChunk$ichunk.gsd --ichunk=$ichunk --tchunk=$tchunk --deltaE=$deltaE\""
	python $exeDIR/BisectChunk.py --user="trajChunk$ichunk.gsd --ichunk=$ichunk --tchunk=$tchunk --deltaE=$deltaE"

	#Now that trajChunk$ichunk has been fully analyzed, we can delete it
	#Uncomment if we want to delete
	#rm -f trajChunk$ichunk.gsd
	
done
#--------------------------------#
# End of the bisection by chunks #
#--------------------------------#

