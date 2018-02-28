#!/bin/bash

declare -A DR_LIST=( ["10.0"]=0.1 ["2.0"]=0.1 ["0.6"]=0.1 )



for T in 10.0 2.0 0.6
do
    for N in 65
    do
	for ISAM in 0
	do
	    conf=../OUTPUT/T$T/N$N/S$ISAM/thermalized.gsd
	    dr=${DR_LIST[$T]}
	    
	    python ../THERMALIZE/progs/CalculatePairCorrelationFunction.py --user="$conf --dr=$dr"
	    mv gofr.txt FIGURES/gofr_T${T}N${N}i${ISAM}.txt
	    mv gofr.png FIGURES/gofr_T${T}N${N}i${ISAM}.png
	done
    done
done

