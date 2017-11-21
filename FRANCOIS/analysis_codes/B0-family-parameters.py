## short script to priduce families of values of the descriptor functions parameters.
from __future__ import print_function
import gsd.pygsd
import gsd.hoomd
import gsd.fl

import sys
import numpy as np

import module_filenamesHandler 
import module_PBC
import module_descriptors

######
## plots: 
import module_importPlotParams
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True          
plt.rcParams.update(module_importPlotParams.fparams())
plt.ion()   ## optional
######

################################################################################

## one is free to define the family of descriptor(s) in terms of functional forms and their parameters ##

descriptor_family_name= 1 # 'family1'

dr = 0.1
r_cutoff= 3.0 + dr/2.
sigma=0.1

r_min = 0.7
rs = np.arange(r_min, r_cutoff, dr)
M_radial_XY = len(rs)
sigmas = np.ones(M_radial_XY) * sigma
descriptor_family_parameters = np.array([rs, sigmas]).transpose()
header = "family of type #1, purely radial parameters, with r_cutoff="+str(r_cutoff)+"   dr="+str(dr)+"   r_min="+str(r_min)

np.savetxt("B0_fam="+str(descriptor_family_name)+'_type=radial.dat', descriptor_family_parameters, fmt="%1.5f %1.5f", header=header)


################################################################################
#descriptor_family_name= 2 # 'family1'

#descriptor_family_parameters = np.array([xis, zetas, lambdas]).transpose()

#np.savetxt("B0_fam="+str(descriptor_family_name)+'_type=radial.dat', descriptor_family_parameters, fmt="%.3f %.3f")


