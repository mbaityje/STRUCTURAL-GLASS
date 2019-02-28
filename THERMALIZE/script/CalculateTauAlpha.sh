#!/bin/bash

# Estimate tau_alpha crudely, and write it in the data folder

thermostat=${thermostat:-'NVE'}


tabella=../data/tau_alpha.txt
echo "T N tau" > $tabella

for T in 5.0 2.0 1.0 0.8 0.7 0.6 0.49 0.46
do
	for N in 65
	do
		Fkt_file="../../OUTPUT/T$T/N$N/shift/Fkt_aftergap_shift_${thermostat}.txt"

		if ! [ -a $Fkt_file ]; then continue; fi

		tau=`python -c "import sys,numpy as np; \
		a=np.loadtxt(sys.argv[1]); \
		i=np.where(a[:,1]<np.exp(-1))[0][0];\
		print(a[i][0])" $Fkt_file`
		echo "$T $N $tau" >> $tabella
	done
done
