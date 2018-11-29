#!/bin/bash

#
# This script calculates the average self-intermediate scattering function.
#
nskip=2 

rootDIR=../..
outDIR=$rootDIR/OUTPUT

TLIST=${1:-"10.0 5.0 2.0 1.0 0.8 0.7 0.6 0.55 0.52 0.49 0.466 0.45 0.44 0.43"}
NLIST=${2:-"1080 65"}

function average {
	awk -vnskip=$nskip '\
							BEGIN{print "#1)t 2-3)Fkt";}
							(NR>nskip){}{if(NR==FNR){t[FNR-nskip]=$1}\
							sum[FNR-nskip]+=$2;\
							sum2[FNR-nskip]+=$2*$2; n[FNR-nskip]++;}
							END\
							{for(i=1;i<=FNR-nskip;i++)\
							{mean=sum[i]/n[i]; mean2=sum2[i]/n[i]; err=sqrt((mean2-mean*mean)/(n[i]-1));\
							print t[i],mean,err;}}' $input
}


for T in $(echo $TLIST)
do
    for N in $(echo $NLIST)
    do
	for tipo in ifr0 aftergap
	do
	    if [ $N == 1080 ];
	    then
		if [ -z $pot_mode ]; then pot_mode=xplor; fi
		input=$outDIR/T$T/N$N/S[0-9]/Fkt_${tipo}_${pot_mode}.txt
		output=$outDIR/T$T/N$N/Fkt_${tipo}_${pot_mode}.txt
	    elif [ $N == 65 ]
	    then
		 if [ -z $pot_mode ]; then pot_mode=shift; fi
		 input=$outDIR/T$T/N$N/$pot_mode/S[0-9]/Fkt_${tipo}_${pot_mode}.txt
		 output=$outDIR/T$T/N$N/$pot_mode/Fkt_${tipo}_${pot_mode}.txt
	    fi
	    average $input > $output
	done
    done
done
