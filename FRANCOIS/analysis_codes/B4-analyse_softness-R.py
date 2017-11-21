### B4
from __future__ import print_function

import gsd.pygsd
import gsd.hoomd
import gsd.fl
import numpy as np
import sys

import module_filenamesHandler 

######
## plots: 
import module_importPlotParams
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True          
plt.rcParams.update(module_importPlotParams.fparams())
plt.ion()   ## optional
######


print('usage:')
print('run  '+sys.argv[0]+'   [filename]   [kind]')

filename = sys.argv[1]
kind     = sys.argv[2]
#    kind = module_filenamesHandler.filename_parser(filename[:-4], 'kind')

if "softness-R" in filename :   # read the decorated file produced before # 

    p_threshold_low = float(module_filenamesHandler.filename_parser(filename[:-4], 'pLow'))
#    print('reading that the file has a p_threshold_low=', p_threshold_low, ' ...')
#    if p_threshold_low ==0 : 
#        p_threshold_low= 1e-7
#        print("correcting its threshold to :",p_threshold_low)

#    Natoms = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'N')) )
#    selectedAtoms = np.arange(Natoms, dtype=int)
#    if 'tagBonusNum' in filename :
#       Natoms = int(module_filenamesHandler.filename_parser(filename[:-4], 'tagBonusNum')) 
#       selectedAtoms = selectedAtoms[:Natoms]


    decoratedFlowRead = gsd.fl.GSDFile(name=filename, mode='rb')   # TODO : use with as instead of close()
    
    atomTypes = decoratedFlowRead.read_chunk(frame=0, name='atomTypes')
    boxParams   = decoratedFlowRead.read_chunk(frame=0, name='configuration.box')   
    Lx = boxParams[0]
    Ninitial = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'N')))
    expectedDensity = float(module_filenamesHandler.filename_parser(filename[:-4], 'rho'))
    selectedAtoms, Natoms, expectedDensity = module_filenamesHandler.selectAtoms_function(atomTypes, kind, Lx, expectedDensity, 0)
    selectedAtoms = np.arange(Ninitial, dtype=int)[selectedAtoms]
#    selectedAtoms = np.arange(Natoms, dtype=int)
#    if 'tagBonusNum' in filename :
#       Natoms = int(module_filenamesHandler.filename_parser(filename[:-4], 'tagBonusNum')) 
#       selectedAtoms = selectedAtoms[:Natoms]
    atomTypes = atomTypes[selectedAtoms]
#    atomIds   = decoratedFlowRead.read_chunk(frame=0, name='atomIds')[selectedAtoms]


    unique_times_number = decoratedFlowRead.nframes
    unique_times    = np.zeros(unique_times_number, dtype=int)
    phop           = np.zeros((unique_times_number, Natoms))
    softness        = np.zeros((unique_times_number, Natoms))
    for tc_index in range(unique_times_number):
        unique_times[tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='tc')[0]
        phop        [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='phop')[selectedAtoms]    #  resultsArray[:,tc_index,1]
        softness    [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='softness')[selectedAtoms]
    decoratedFlowRead.close()
    
    
################################################################################
    
    
    ### plots ####
    mask = phop>0
    decent_phop = phop[mask]
#    plt.plot()
    
    
    ###### distribution of p_hop_value ######
    plt.figure(40,[20,6])
    base = 1.1
    Log_bins = np.array([p_threshold_low*base**i for i in range(-5, int(np.log(np.max(decent_phop)/p_threshold_low)/np.log(base)), 1) ])
    heights, trash = np.histogram(decent_phop, bins=Log_bins)
    heights = heights / np.diff(Log_bins)    ## adjust density by histo bin width
    plt.loglog(Log_bins[:-1], heights, lw=3)
    plt.xlabel(r'$p_{hop}$')
    plt.ylabel(r'$N(p_{hop} | p_{hop}>'+str(p_threshold_low)+')$')

   
    outName=filename[:-4] + "distro-phop-loglog_phop-Softness"

    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())    
    
    
    plt.figure(41,[20,6])
    heights, trash = np.histogram(softness[mask], bins=100)
#    heights = heights / np.diff(Log_bins)    ## adjust density by histo bin width
    plt.plot(trash[:-1], heights, lw=3)
#    plt.xlabel(r'$p_{hop}$')
#    plt.ylabel(r'$N(p_{hop} | p_{hop}>'+str(p_threshold_low)+')$')
    
    
    scatterplot : 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
