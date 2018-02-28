#!/usr/bin/python

################################################################
#
# DESCRIPTION
# This example reads a configuration from a .gsd file, and
# calculates its energy with several different potentials.
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
hoomd.context.initialize() #Reads arguments from sys.argv and uses them automatically
hoomd.option.set_notice_level(0) #0: no notices - 10: lots of stuff
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd

################################################################
#
# FUNCTIONS THAT WOULD GO IN A SEPARATE MODULE
# 
################################################################

#Returns system's energy if the analyzer is set
def PotEn():
	modeT.set_params(dt=1e-18)
	hoomd.run(2)
	U=analyzer.query('potential_energy')
	return U


################################################################
#
# READ ARGUMENTS
# Read command line arguments
#
################################################################
print("\n\n\nREAD ARGUMENTS\n")

#The we create a parser for the User files 
parser = argparse.ArgumentParser(prog='python ReadAndEvolve.py [--hoomdi-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"',
								 description='The program reads a .gsd configuration and runs it in NVE and NVT',
								 add_help=True)
parser.add_argument('filename', #positional argument
					nargs=1,
					help='name of the .gsd configuration we want to read'
)
args = parser.parse_args(more_arguments)
filename=args.filename[0]
T=0.5
thermostat='NVT'
del parser
print("filename: ",filename)
print("T=",T)
print("Thermostat: ", thermostat)
system = hoomd.init.read_gsd(filename=filename)
Natoms = len(system.particles)
print("Natoms = ",Natoms,"\n")

################################################################
# 
# CLASSES FOR THE POTENTIAL
#
################################################################

class LJ():
	""" A class that contains all the information on the potential
	type: KA, KAshort
	"""
	def __init__(self, NeighborsList, type="KAshort", mode="xplor"):
		self.mode=mode
		assert(self.mode=="xplor") #I only implement xplor, cause it's the only one I use. To obtain the other modes you only need to remove pieces and look at tutorials T11, T12.
		if type=="KA": #Kob-Andersen potential
			self.setKA(NeighborsList)
		elif type=="KAshort": #Kob-Andersen for small systems, as proposed by Heuer
			self.setKA(NeighborsList, r_cutoff=1.4, r_buff=0.0)
		else:
			print("type Non implementato")
			raise SystemExit
 
	def setKA(self,NeighborsList, r_on_cutoff=1.2, r_cutoff=2.5, r_buff=None):
		self.eps_AA=1
		self.eps_AB=1.5
		self.eps_BB=0.5
		self.sig_AA=1
		self.sig_AB=0.8
		self.sig_BB=0.88
		self.r_on_cutoff=r_on_cutoff
		self.r_cutoff=r_cutoff
		self.ron_AA=self.r_on_cutoff*self.sig_AA
		self.ron_AB=self.r_on_cutoff*self.sig_AB
		self.ron_BB=self.r_on_cutoff*self.sig_BB
		self.rcut_AA=self.r_cutoff*self.sig_AA
		self.rcut_AB=self.r_cutoff*self.sig_AB
		self.rcut_BB=self.r_cutoff*self.sig_BB
		if r_buff!=None:
			self.r_buff=r_buff
			NeighborsList.set_params(r_buff=self.r_buff)
		self.myLjPair = md.pair.lj(r_cut=self.r_cutoff, nlist=NeighborsList)
		self.myLjPair.pair_coeff.set('A', 'A', epsilon=self.eps_AA, sigma=self.sig_AA, r_cut=self.rcut_AA, r_on=self.ron_AA)
		self.myLjPair.pair_coeff.set('A', 'B', epsilon=self.eps_AB, sigma=self.sig_AB, r_cut=self.rcut_AB, r_on=self.ron_AB)
		self.myLjPair.pair_coeff.set('B', 'B', epsilon=self.eps_BB, sigma=self.sig_BB, r_cut=self.rcut_BB, r_on=self.ron_BB)
		self.myLjPair.set_params(mode=self.mode)

	def GetLJpair(self):
		return self.myLjPair




class JammingSphere():
	""" A class that contains all the information on the potential
	type: Harmonic, Hertzian
	"""
	def __init__(self, NeighborsList, type="Harmonic", width=10000):
		self.width=width
		self.table=md.pair.table(width=self.width, nlist=NeighborsList, name=type)
		if type=="Harmonic":
			self.alpha=2.0
		elif type=="Hertzian":
			self.alpha=2.5
		else:
			print("type Non implementato")
			raise SystemExit
		self.setCoeff(NeighborsList)

	def JamSphere(self,r,rmin,rmax, radi, radj, alpha=2, V0=1):
		sumrad=radi+radj
		V= (self.V0/self.alpha)*(1-r/sumrad)**self.alpha
		F=-(self.V0/sumrad)    *(1-r/sumrad)**(self.alpha-1.0)
		return (V,F)

	def setCoeff(self,NeighborsList):
		self.mode='None'
		self.V0=1
		self.radA=0.5303030303030303 #Chose rA and rB such that in an 80:20 mixture the average sphere diameter is 1, and rA/rB=1.4
		self.radB=0.3787878787878788
		self.table.pair_coeff.set('A','A',func=self.JamSphere, rmin=0, rmax=2*self.radA,		coeff=dict(radi=self.radA, radj=self.radA, alpha=self.alpha, V0=self.V0))
		self.table.pair_coeff.set('A','B',func=self.JamSphere, rmin=0, rmax=self.radB+self.radA,coeff=dict(radi=self.radA, radj=self.radB, alpha=self.alpha, V0=self.V0))
		self.table.pair_coeff.set('B','B',func=self.JamSphere, rmin=0, rmax=2*self.radB,		coeff=dict(radi=self.radB, radj=self.radB, alpha=self.alpha, V0=self.V0))
		# self.table.pair_coeff.set('A','A',func=self.JamSphere, rmin=0, rmax=2*self.radA, coeff=dict())
		# self.table.pair_coeff.set('A','B',func=self.JamSphere, rmin=0, rmax=self.radB+self.radA, coeff=dict())
		# self.table.pair_coeff.set('B','B',func=self.JamSphere, rmin=0, rmax=2*self.radB, coeff=dict())

	def GetPair(self):
		return self.table

	def GetWidth(self):
		""" 
		How many points there are in the table that defines the potential. More points=more precise energy
		"""
		return self.width


#Returns a hoomd.md.pair with Kob-Anderson Binary Lennard-Jones parameters
def KApotential(NeighborsList, r_cutoff=2.5): 
	eps_AA=1 
	eps_AB=1.5
	eps_BB=0.5
	sig_AA=1
	sig_AB=0.8
	sig_BB=0.88
	r_on_cutoff=1.2
	myLjPair = md.pair.lj(r_cut=r_cutoff, nlist=NeighborsList)
	myLjPair.pair_coeff.set('A', 'A', epsilon=eps_AA, sigma=sig_AA, r_cut=r_cutoff*sig_AA, r_on=r_on_cutoff*sig_AA)
	myLjPair.pair_coeff.set('A', 'B', epsilon=eps_AB, sigma=sig_AB, r_cut=r_cutoff*sig_AB, r_on=r_on_cutoff*sig_AB)
	myLjPair.pair_coeff.set('B', 'B', epsilon=eps_BB, sigma=sig_BB, r_cut=r_cutoff*sig_BB, r_on=r_on_cutoff*sig_BB)
	myLjPair.set_params(mode="xplor")
	return myLjPair

#Returns a hoomd.md.pair with Kob-Anderson Binary Lennard-Jones parameters
def KApotentialShort(NeighborsList):
	eps_AA=1
	eps_AB=1.5
	eps_BB=0.5
	sig_AA=1
	sig_AB=0.8
	sig_BB=0.88
	r_on_cutoff=1.2
	r_cutoff=1.4
	NeighborsList.set_params(r_buff=0.0)
	myLjPair = md.pair.lj(r_cut=r_cutoff, nlist=NeighborsList)
	myLjPair.pair_coeff.set('A', 'A', epsilon=eps_AA, sigma=sig_AA, r_cut=r_cutoff*sig_AA, r_on=r_on_cutoff*sig_AA)
	myLjPair.pair_coeff.set('A', 'B', epsilon=eps_AB, sigma=sig_AB, r_cut=r_cutoff*sig_AB, r_on=r_on_cutoff*sig_AB)
	myLjPair.pair_coeff.set('B', 'B', epsilon=eps_BB, sigma=sig_BB, r_cut=r_cutoff*sig_BB, r_on=r_on_cutoff*sig_BB)
	myLjPair.set_params(mode="xplor")
	return myLjPair


################################################################
# 
# SET UP INITIAL POTENTIAL
#
################################################################
NeighborsList = md.nlist.cell()
analyzer_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']
analyzer = hoomd.analyze.log(filename=None, quantities=analyzer_quantities, period=1)
modeT=md.integrate.mode_standard(dt=1e-18)
integrator=md.integrate.nve(group=hoomd.group.all())


print("KOB-ANDERSEN SHORT POTENTIAL")
KAshort=KApotentialShort(NeighborsList)
print("E=",PotEn())
KAshort.disable()

print("KOB-ANDERSEN SHORT POTENTIAL WITH CLASS")
potKAshort=LJ(NeighborsList,type="KAshort") #myLJpair is now an attribute of potential. To call it: potential.GetLJpair()
print("E=",PotEn())
potKAshort.GetLJpair().disable()

if Natoms>500:
	print("KOB-ANDERSEN POTENTIAL")
	KA=KApotential(NeighborsList)
	print("E=",PotEn())
	KA.disable()

	print("KOB-ANDERSEN POTENTIAL WITH CLASS")
	potKA=LJ(NeighborsList,type="KA") #myLJpair is now an attribute of potential. To call it: potential.GetLJpair()
	print("E=",PotEn())
	potKA.GetLJpair().disable()


print("HARMONIC SPHERES")
def JamSphere(r,rmin,rmax,radi, radj, alpha=2.0, V0=1):
	sumrad=radi+radj
	V= (V0/alpha)*(1-r/sumrad)**alpha
	F=-(V0/sumrad)*(1-r/sumrad)**(alpha-1.0)
	return (V,F)

harm=md.pair.table(width=1000, nlist=NeighborsList, name='harmonicSpheres')
radA=0.5303030303030303 #Chose rA and rB such that the average sphere diameter is 1, and rA/rB=1.4
radB=0.3787878787878788
alpha=2.0
V0=1
harm.pair_coeff.set('A','A',func=JamSphere, rmin=0, rmax=2*radA,    coeff=dict(radi=radA, radj=radA, alpha=alpha, V0=V0))
harm.pair_coeff.set('A','B',func=JamSphere, rmin=0, rmax=radB+radA, coeff=dict(radi=radA, radj=radB, alpha=alpha, V0=V0))
harm.pair_coeff.set('B','B',func=JamSphere, rmin=0, rmax=2*radB,    coeff=dict(radi=radB, radj=radB, alpha=alpha, V0=V0))
print("E=",PotEn())
harm.disable()



print("HARMONIC SPHERES THROUGH CLASS")
potHarm=JammingSphere(NeighborsList,type="Harmonic")
print("E=",PotEn())
potHarm.GetPair().disable()


