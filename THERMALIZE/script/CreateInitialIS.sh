#!/bin/bash

rootDIR=../..

exeDIR=$PWD/../progs
cd $exeDIR
iniDIR=$rootDIR/OUTPUT/INITIAL-STATES/N65/
mkdir -p $iniDIR

cd $iniDIR
python $exeDIR/CreateInitialIS_N65.py




