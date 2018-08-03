#!/usr/bin/env python
################################################################
#
#
# DESCRIPTION
# This example creates many high-temperature configurations and minimizes them. The idea is to see the distribution of the inherent-structure energies, to see if we can define a threshold energy.
#
# To display help:
# python T24-ThresholdEnergy.py --user="-h"
#
# To launch a simulation:
# python T24-ThresholdEnergy.py
#
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md #md is the subsection strictly about the properties of the system
import sys #this has some system tools like sys.exit() that can be useful
import argparse #for processing arguments
import numpy as np #Handles some mathematical operations
from math import sqrt #Just wanna square root
import gsd.pygsd #gsd is the database type for the configurations
import gsd.hoomd #
from lib import module_potentials as pot


def Minimize(snap):
	system.restore_snapshot(snap)
	fire.cpp_integrator.reset()
	while not(fire.has_converged()):
		hoomd.run(100)
	eIS=analyzer.query('potential_energy')
	return eIS


################################################################
#
# SET UP THE SIMULATION CONTEXT
# 
################################################################
print("\n\n\nSET UP THE SIMULATION CONTEXT\n")
hoomd.context.initialize() #Reads arguments from sys.argv and uses them automatically
hoomd.option.set_notice_level(0) #0: no notices - 10: lots of stuff
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd

#Read parameters
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename', help='name of the trajectory file .gsd')
parser.add_argument('-n','--nconfs', type=int, required=False, default=10, help='number of IS to generate')
parser.add_argument('--dt', type=float, required=False, default=0.0025, help='time step')
args = parser.parse_args(more_arguments)




#Load configuration
system = hoomd.init.read_gsd(filename=args.filename)
Natoms=len(system.particles)
#Set Potential
NeighborsList = md.nlist.cell()
mypot="KAshort" if Natoms<500 else "KA"
potential=pot.LJ(NeighborsList,type=mypot)
analyzer = hoomd.analyze.log(filename=None, quantities=['temperature','potential_energy', 'kinetic_energy', 'momentum'], period=int(1./args.dt), header_prefix = '#', overwrite=True, phase=0)


eislist=[]
eTlist=[]
for i in range(args.nconfs):
	modeT=md.integrate.mode_standard(dt=args.dt)
	if i==0: integratorNVT = md.integrate.nvt(group=hoomd.group.all(), kT=10.0, tau=0.1)
	else: integratorNVT.enable()
	integratorNVT.randomize_velocities(seed=int( (i*1664525+1013904223))%4294967295)
	hoomd.run(2000,quiet=True)
	eT=analyzer.query('potential_energy')
	integratorNVT.disable()

	fire=hoomd.md.integrate.mode_minimize_fire(dt=0.0025, alpha_start=0.99, ftol=1e-5, Etol=1e-10, wtol=1e-5)
	if i!=0: integratorFIRE.enable()
	else: integratorFIRE = md.integrate.nve(group=hoomd.group.all())
	snap=system.take_snapshot()
	eis=Minimize(snap)
	integratorFIRE.disable()
	print('i: ',i,'eT= ',eT,' eis= ',eis)
	eTlist.append(eT)
	eislist.append(eis)



#Visualize distribution

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

sns.distplot(eislist)
sns.distplot(eTlist)
plt.show()
sns.distplot(eislist)
plt.show()
df = pd.DataFrame(np.array([eislist,eTlist]).transpose(), columns=["x", "y"])
sns.jointplot(x='x',y='y',data=df, kind="kde");
plt.show()
