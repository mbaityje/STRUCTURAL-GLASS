#!/bin/bash
#PBS -V
#PBS -q longgpu
#PBS -M francois.landes@gmail.com
#PBS -m bae
#PBS -N n64-T35-long
#PBS -l walltime=168:00:00
#PBS -e logs/$ro.$PBS_JOBID.e.txt
#PBS -o logs/$ro.$PBS_JOBID.o.txt

###########################################################

###############################
### parameters ################
ppn=1
days=7

# T45-27-glued-disappear-peak
## 1e6 at 0.45 then 1e7 at NVE 
# T45-scratch-Heat045
## cool-dup-045to043


rootname=$ro
#rootname="KA_rho=12e-1_N=1000_T=45e-2_tag=196-10"
#rootname="KA_rho=12e-1_N=1000_T=45e-2_tag=196-10_dt=0.0008"

echo "need arguemnts: $ro  "
 
### T=0.41
# qsub -N T41ultraLongNVT    -v ro=KA_rho\=12e-1_N\=8000_T\=41e-2_tag\=900-12345678  sub.sh
# tstepsNVT=20e8
# recPeriod=1e8
# time = 336h

### T=0.45
# qsub -v ro=KA_rho\=12e-1_N\=1000_T\=45e-2_tag\=196-10  sub.sh

######################
### T=0.43
# qsub -v ro=KA_rho\=12e-1_N\=1000_T\=43e-2_tag\=180-2   sub.sh

# qsub -v ro=KA_rho\=12e-1_N\=8000_T\=43e-2_tag\=scr1   sub.sh
## remember to change T=5 to 1 to ## and then only run glueing_flag=glueNVT
# NVT-43-scr1-initialNVT
# NVT-43-scr1-init (==frame 0)

# qsub              -v ro=KA_rho\=12e-1_N\=8000_T\=43e-2_tag\=180-12345678  sub.sh 
# qsub -N glueMore  -v ro=KA_rho\=12e-1_N\=8000_T\=43e-2_tag\=180-12345678  sub.sh

# stitched-43-1-8-NVTglue
# qsub -v ro=KA_rho\=12e-1_N\=8000_T\=43e-2_tag\=180-12345678-Fr-100   sub.sh
# NVE-43-Fr999-glu
# NVE-43-Fr100-glu

#NVE-43-Fr999-longerNVE

restartFileName=$rootname\_type=restart-initial.gsd
restartFileName=$rootname\_type=restart.gsd

#tstepsNVT=20e8 #2e8
#tol_NVT_conv=0
#tstepsNVE=0
#glueing_flag=glueNVT
#recPeriod=100000000 #2e7
#analyzer_period=1e6
#tauT=0.01


tstepsNVT=20e7  # in 64000 , 20e7 should be 4.6 days 
tol_NVT_conv=-1 # 0.001 # or -1 (shift) or -2 (xplor)
tstepsNVE=1e6 # 1e6 # 3000500 #1e8  ## 4e6 = 40 min in N=8000 . 
glueing_flag=0
recPeriod=4000
analyzer_period=4000
tauT=1000.0
dt=0.0025 # $dtArg # 0.001
thermostat=NVT
## end of parameters section ##
###############################



##########################################################
## directory from which this job is submitted
cd $PBS_O_WORKDIR   

# create directories on the node local HDD 
mkdir /scratch/flandes
TMP=/scratch/flandes/$USER.$PBS_JOBID
mkdir $TMP
cp $PBS_O_WORKDIR/$rootname*   $TMP/.
cp $PBS_O_WORKDIR/*.py         $TMP/.
cd $TMP



logFileName=$rootname\_type=XPU-times.log
touch $logFileName
command_ran_name=$rootname\_type=command_ran.sh
echo " " >  $command_ran_name
#restartFileName="$rootname\_type=restart.gsd"
export HOOMD_WALLTIME_STOP=$(( `date +%s` + $days * 24 * 3600 - 2*3600 )) ## 4 days minus 2 hours
echo $HOOMD_WALLTIME_STOP
#module list 
#module load pythonXXX


#device=0
# XXX DEBUG OPTIONS: mpirun -np 2 cuda-memcheck hoomd myscript.py —gpu_error_checking
## --mode=gpu  
## --notice-level=0

echo "time  mpirun -n $ppn   python   MD1-run_NVT_T_NVE.py   $restartFileName  $tstepsNVT  $tol_NVT_conv  $tstepsNVE   $glueing_flag   $recPeriod   $analyzer_period  $tauT  $dt $thermostat  >> $logFileName " > $command_ran_name

echo " " >> $command_ran_name
chmod +x $command_ran_name

$(./$command_ran_name) & pid=$!



#device=1
##seed=$(($seed_MaxAttempts+$device))
##dprime=$(echo "$multiplier * $d" | bc -l )
#echo "" > $command_ran_name
#chmod +x $command_ran_name
#$(./$command_ran_name) & pid2=$!
#wait $pid2 


## wating for both indep. runs to finish:
wait $pid 
#wait $pid2

######## end of simulation here ##########

## take back data.. be very careful to test if you change it.
mv $rootname\_type=* $PBS_O_WORKDIR/.
mv $rootname\_dt=*\_type=* $PBS_O_WORKDIR/.
cd ..
mv $TMP $PBS_O_WORKDIR/logs/.
cp  $PBS_O_WORKDIR/$rootname\_type=restart.gsd   $PBS_O_WORKDIR/logs/$USER.$PBS_JOBID/.
cp  $PBS_O_WORKDIR/$rootname\_type=notes.txt     $PBS_O_WORKDIR/logs/$USER.$PBS_JOBID/.

exit 



################################################################################
################################################################################
#module load intel-13.1
### netcdf library- 
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/pkg/netcdf-cxx/4.2-icc/lib:/opt/pkg/netcdf/4.3.0-icc/lib
###hdf library
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/pkg/hdf5/1.8.10-patch1-icc/lib
#../xidentifyHopsv2 Inputs > out.identifiyhops

#date > out.construct
#../xconstructTrainingSetHopsv2_Ni Inputs >> out.construct

#date > out.trainsvm
#../xtrainSVM Inputs >> out.trainsvm

#../xassignSoftnessField Inputs > out.assignsoft



### T=5 and 1 - for annealing ##
#tstepsNVT=1e7
#tol_NVT_conv=0
#tstepsNVE=0
#glueing_flag=0
#recPeriod=1e3
#analyzer_period=1e4
#tauT=0.01


