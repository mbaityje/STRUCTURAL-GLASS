#!/bin/bash
#SBATCH --ntasks=1
#SBATCH -p gpu # partition (queue)
#SBATCH --gres=gpu:1

# 
# Bisects a chunk of trajectory with Heuer's IS bisection scheme.
# 


if [ `hostname` == "PennPuter" ] || [ `hostname` == "tango" ]; then SYSTEM="PennPuter";
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
scriptDIR=$thermDIR/script
workDIR=$rootDIR/OUTPUT
utilDIR=$rootDIR/UTILITIES



#
#Command line input. Default values are for testing purposes. 
#
T=${T:-$1}
dt=${dt:-$2}
deltaE=${deltaE:-$3}
ttot=${ttot:-$4}
doridge=${doridge:-$5}
doridge=${doridge:-0} #default value for doridge


#
#Some checks to make sure that the input is good
#
#Number of arguments
if [ $SYSTEM == "PennPuter" ] && [ $# -ne 4 ] && [ $# -ne 5 ]; then
	echo "Wrong number of parameters ($#)"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doridge>]"
	echo "T: temperature"
	echo "dt: MD integration step"
	echo "deltaE: energy difference to state that two neighboring IS are different"
	echo "ttot: total length of the trajectory"
	echo "doridge: 0 (default): Calculate the trajectory of the IS. 1: Calculate both IS and TS."
	exit
fi
#T must be a non negative number
if ! [[ $T =~ ^[0-9]*\.?[0-9]*$ ]]; then
	echo "T (=$T) must be a non-negative number"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doridge>]"
	exit
fi
#dt must be 0<dt<0.01
if ! [[ $dt =~ ^[0-9]*\.?[0-9]*$ ]]; then
	echo "dt (=$deltaE) must be a non-negative number"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doridge>]"
	exit
elif [ 1 -eq $(echo "$dt > 0.01" | bc -lq)  ]; then
	echo "dt = $dt >= 0.01, is very big, you will likely have particles falling out of the box"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doridge>]"
	exit
fi
#deltaE must be a non negative number
if ! [[ $deltaE =~ ^[0-9]*\.?[0-9]*$ ]]; then
	echo "deltaE (=$deltaE) must be a non-negative integer"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doridge>]"
	exit
fi
#ttot must be a positive integer
if ! [[ $ttot =~ ^[0-9]+$ ]]; then
	echo "ttot (=$ttot) must be a positive integer"
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doridge>]"
	exit
fi
#doridge must be 0 or 1
if ! [ $doridge -eq 0 -o $doridge -eq 1 ]; then
	echo "doridge ($doridge) must be 0 or 1."
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doridge>]"
	exit
fi

#doridge must be 0 or 1
if [ $doridge -eq 0 ]; then
	doridgeflag=''
elif [ $doridge -eq 1 ]; then
	doridgeflag='--doridge'
else
	echo "doridge ($doridge) must be 0 or 1."
	echo "Launch as: $0 <T> <dt> <deltaE> <ttot> [<doridge>]"
	exit
fi


#Hard-coded parameters
readonly tauT=0.1
readonly tchunk=`echo 10^5|bc` #`echo 10^5|bc`
nchunks=`echo $ttot/$tchunk|bc`
let nchunksm1=$nchunks-1

#--------------------------------------------#
# Here starts the actual bisection by chunks #
#--------------------------------------------#

# Filenames
filename=thermalized.gsd
elistFILE=elistIS.txt

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
echo "elistFILE: $elistFILE"
echo "firstchunk: = $firstChunk"
echo "nchunks: = $nchunks"
echo "tlast: = $tlast"
for ichunk in $(seq $firstChunk $nchunksm1)
do
	echo "+++ ichunk = $ichunk +++"
	STARTTIME=$(date +%s)

	#First generate the thermal trajectory of the cunk
	python $exeDIR/CreateChunk.py --user="$filename --ichunk=$ichunk --tchunk=$tchunk --dt=$dt --temperature=$T --tauT=$tauT"

	#Now that we have the new chunk, we can delete the previous one
	#Uncomment the next line if we don't want to keep the thermal trajectory
	#if [ $ichunk -ge 1 ]; then rm -f restartChunk`expr $ichunk - 1`.gsd; fi

	#Perform IS bisection on the chunk
	if [ $ichunk == 0 ]; then skiprows=0;
	elif ! [ -f $elistFILE ]; then echo "$elistFILE does not exist. Should exist since ichunk=$ichunk"; exit;
	else skiprows=`wc -l $elistFILE | awk '{print $1-2}'`; fi
	echo "skiprows: $skiprows"

	echo "python $exeDIR/BisectChunkWithRidge.py --user=\"trajChunk$ichunk.gsd --ichunk=$ichunk --tchunk=$tchunk --deltaE=$deltaE $doridgeflag --skiprows=$skiprows\""
	python $exeDIR/BisectChunkWithRidge.py --user="trajChunk$ichunk.gsd --ichunk=$ichunk --tchunk=$tchunk --deltaE=$deltaE $doridgeflag --skiprows=$skiprows"

	#Now that trajChunk$ichunk has been fully analyzed, we can delete it
	#Uncomment if we want to delete
	#rm -f trajChunk$ichunk.gsd

	ENDTIME=$(date +%s)
	ELAPSEDTIME=`echo "($ENDTIME - $STARTTIME)/60"|bc -l`
	echo "ichunk: $ichunk  time: $ELAPSEDTIME min >> bisect.profile"
done

#--------------------------------#
# End of the bisection by chunks #
#--------------------------------#

