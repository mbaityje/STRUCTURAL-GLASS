#!/bin/bash

#
# This script calculates the average self-intermediate scattering function.
#
nskip=2 

rootDIR=../..
outDIR=$rootDIR/OUTPUT


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

for T in 10.0 2.0 0.6
do
    for N in 65
    do
	for tipo in ifr0 aftergap
	do
	    input=$outDIR/T$T/N$N/S*/Fkt_${tipo}.txt
	    output=$outDIR/T$T/N$N/Fkt_${tipo}.txt
	    average $input > $output
	done
    done
done
