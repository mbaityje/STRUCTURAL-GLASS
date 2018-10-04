#!/usr/bin/env python
################################################################
#
# DESCRIPTION
# This example reads a configuration from a .gsd file, and
# let's it evolve with NVT (first) dynamics, comparing the difference
# when the tau of the thermostat is changed..
#
# To display help:
# python T31-CompareTau.py --user="-h"
#
# To launch a simulation:
# python T31-CompareTau.py --user="filename -t runsteps"
#
# For example:
# python T31-CompareTau.py --user="sample-states/KA_rho=12e-1_N=1e3_T=45e-2_tag=1_type=restart.gsd -e 10 -t 1000 --thermostat NVT"
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md #md is the subsection strictly about the properties of the system
import sys #this has some system tools like sys.exit() that can be useful
import argparse #for processing arguments
import numpy as np #Handles some mathematical operations
from lib import module_potentials as pot
import pandas as pd


hoomd.context.initialize() #Reads arguments from sys.argv and uses them automatically
hoomd.option.set_notice_level(0) #0: no notices - 10: lots of stuff
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd

parser = argparse.ArgumentParser(prog='python ReadAndEvolve.py [--hoomdi-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
parser.add_argument('filename', help='name of the .gsd configuration we want to read')
parser.add_argument('-t','--runsteps', type=int, required=False, default=1000, help='number of NVT steps (default: 1000)')
parser.add_argument('-p','--pot_type', required=False, default='KAshort', help='Type of potential (default: KAshort)')
parser.add_argument('-T','--temperature', type=float, required=False, default=0.45, help='temperature (default: 2.0)')
args = parser.parse_args(more_arguments)

system = hoomd.init.read_gsd(filename=args.filename)
inistep=hoomd.get_step()
ini_snap= system.take_snapshot()
potential=pot.LJ(md.nlist.cell(), type=args.pot_type)
quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum']


md.integrate.mode_standard(dt=0.0025) ##allows Nose-Hoover in NVT, for instance
md.update.zero_momentum(phase=0, period=400)

listatauT=[0.005, 0.01, 0.1, 1.0, 2.0]
lognames=["test-output/test-tau"+str(tauT)+".txt" for tauT in listatauT]


from matplotlib import pyplot as plt

for i in range(len(listatauT)):
	system.restore_snapshot(ini_snap)
	logname=lognames[i]
	analyzer = hoomd.analyze.log(filename=logname, quantities=quantities, period=100, overwrite=True, phase=0)

	integrator_nvt = md.integrate.nvt(group=hoomd.group.all(), kT=args.temperature, tau=listatauT[i])
	hoomd.run(args.runsteps, quiet=False) 
	integrator_nvt.disable()
	analyzer.disable()

	df=pd.read_csv(lognames[i], delim_whitespace=True)
	plt.plot( df['timestep']-inistep, df['potential_energy'], label='tau = '+str(listatauT[i]) )
	inistep=hoomd.get_step()

plt.legend()
plt.show()


