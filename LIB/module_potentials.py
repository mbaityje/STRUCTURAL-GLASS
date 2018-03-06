#!/usr/bin/python
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md
import numpy as np
import os.path
if os.path.isfile('module_measurements.py'):
	import module_measurements as med #python2 import style
else:
	import lib.module_measurements as med #python3 import style

class LJ:
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

	#XPLOR smoothing function
	def S(self, r, r_on, r_cut):
		assert(r>0)
		if r<r_on:
			return 1
		elif r>r_cut:
			return 0
		else:
			return self.Stilde(r, r_on, r_cut)
	#Stilde is the non-trivial part of S
	def Stilde(self, r, r_on, r_cut):
		assert(r>r_on)
		assert(r<r_cut)
		rc2=r_cut*r_cut
		ro2=r_on*r_on
		r2=r*r
		term1=rc2 - r2
		term2=rc2 + 2*r2 - 3*ro2
		term3=rc2 - ro2
		return (term1*term1*term2)/(term3*term3*term3)

	#For the derivative of S, I only consider the non-trivial part
	def StildePrime(self, r, r_on, r_cut):
		assert(r>r_on)
		assert(r<r_cut)
		rc2=r_cut*r_cut
		ro2=r_on*r_on
		r2=r*r
		term1=rc2 - r2
		term2=ro2 - r2
		term3=rc2 - ro2
		return 12*r*term1*term2/(term3*term3*term3)
	#LJ pair pure potential
	def Vpair(self, r, eps, sigma, r_cut):
		assert(r>0)
		if r>=r_cut:
			return 0
		else:
			sigma_on_r=sigma/r
			temp=sigma_on_r*sigma_on_r
			sigma_on_r6=temp*temp*temp
			return 4*eps*(sigma_on_r6*sigma_on_r6 - sigma_on_r6)
	#Derivative of the LJ pair pure potential
	def VpairPrime(self, r, eps, sigma, r_cut):
		if r>=r_cut:
			return 0
		invr=1./r
		invr2=invr*invr
		invr4=invr2*invr2
		invr6=invr4*invr2
		invr8=invr4*invr4
		s2=sigma*sigma
		s6=s2*s2*s2
		return 48*eps*s6*invr8*(0.5 - s6*invr6)


	#LJ potential with the XPLOR smoothing
	def Vij(self, r, eps, sigma, r_on, r_cut):
		if r>=r_cut:
			return 0
		if(r_on<r_cut):
			return self.S(r,r_on,r_cut)*self.Vpair(r,eps, sigma,r_cut)
		elif(r_on>=r_cut):
			return self.Vpair(r,eps,sigma,r_cut)-self.Vpair(r_cut,eps,sigma,r_cut)

	#Derivative of the potential with XPLOR smoothing
	def Vprime(self, r, eps, sigma, r_on, r_cut):
		if r>=r_cut:
			return 0
		elif r<=r_on:
			return self.VpairPrime(r, eps, sigma, r_cut)
		else:
			return self.VpairPrime(r, eps, sigma, r_cut)*self.Stilde(r,r_on,r_cut)+(self.Vpair(r, eps, sigma, r_cut)*self.StildePrime(r, r_on, r_cut))/r

	def Cd(self,snapA, snapB, beta):
		"""
		Calculate Diagonal part of the correlator defined by Szamel
		"""
		Natoms=snapA.particles.N
		assert(Natoms==snapB.particles.N)
		L=snapA.box.Lx #Assume cubic box
		assert(L==snapB.box.Lz)
		# gradVector = np.zeros((Natoms, 3)) 
		invsqrt3=1./np.sqrt(3) #Assume three dimensions
		Kvec=np.array([invsqrt3 for i in range(3)],dtype=np.float64)

		output=0
		for k in range(Natoms):
			posAk=np.array(snapA.particles.position[k], dtype=np.float64)
			posBk=np.array(snapB.particles.position[k], dtype=np.float64)
			typeAk=snapA.particles.typeid[k]
			typeBk=snapB.particles.typeid[k]
			assert(typeAk==typeBk)
			for i in range(Natoms):
				if i==k:
					continue
				posAi=np.array(snapA.particles.position[i], dtype=np.float64)
				posBi=np.array(snapB.particles.position[i], dtype=np.float64)
				typeAi=snapA.particles.typeid[i]
				typeBi=snapB.particles.typeid[i]
				assert(typeAi==typeBi)
				rvecA=med.PeriodicDisplacement(posAk,posAi,L)
				rvecB=med.PeriodicDisplacement(posBk,posBi,L)
				rA=np.linalg.norm(rvecA)
				rB=np.linalg.norm(rvecB)
				if typeAi==typeAk:
					if typeAi==0:
						VpA=self.Vprime(rA, self.eps_AA, self.sig_AA, self.ron_AA, self.rcut_AA) if rA<=self.rcut_AA else 0
						VpB=self.Vprime(rB, self.eps_AA, self.sig_AA, self.ron_AA, self.rcut_AA) if rB<=self.rcut_AA else 0
					else:
						VpA=self.Vprime(rA, self.eps_BB, self.sig_BB, self.ron_BB, self.rcut_BB) if rA<=self.rcut_BB else 0
						VpB=self.Vprime(rB, self.eps_BB, self.sig_BB, self.ron_BB, self.rcut_BB) if rB<=self.rcut_BB else 0
				else:
					VpA=self.Vprime(rA, self.eps_AB, self.sig_AB, self.ron_AB, self.rcut_AB) if rA<=self.rcut_AB else 0
					VpB=self.Vprime(rB, self.eps_AB, self.sig_AB, self.ron_AB, self.rcut_AB) if rB<=self.rcut_AB else 0
				FA=rvecA*VpA
				FB=rvecB*VpB
				output+=FA.dot(Kvec)*FB.dot(Kvec)
		return output*beta/Natoms


	def CalculateGradient(self,snapshot):
		"""
		Calculate INTENSIVE Energy and Gradient of the smoothened LJ potential
		"""
		energia=0
		Natoms=snapshot.particles.N
		L=snapshot.box.Lx #Assume cubic box
		gradVector = np.zeros((Natoms, 3)) #Assume three dimensions
		for i in range(Natoms):
			for j in range(i+1,Natoms):
				posi=np.array(snapshot.particles.position[i], dtype=np.float64)
				posj=np.array(snapshot.particles.position[j], dtype=np.float64)
				typei=snapshot.particles.typeid[i]
				typej=snapshot.particles.typeid[j]
				rvec=med.PeriodicDisplacement(posi,posj,L)
				r=np.sqrt(np.square(rvec).sum())

				#See what group the particles belong to
				if typei==typej:
					if typej==0:
						if r>self.rcut_AA:
							continue
						V=self.Vij(r, self.eps_AA, self.sig_AA, self.ron_AA, self.rcut_AA)
						Vp=self.Vprime(r, self.eps_AA, self.sig_AA, self.ron_AA, self.rcut_AA)
					else:
						if r>self.rcut_BB:
							continue
						Vp=self.Vprime(r, self.eps_BB, self.sig_BB, self.ron_BB, self.rcut_BB)
						V=self.Vij(r, self.eps_BB, self.sig_BB, self.ron_BB, self.rcut_BB)
				else:
					if r>self.rcut_AB:
						continue
					Vp=self.Vprime(r, self.eps_AB, self.sig_AB, self.ron_AB, self.rcut_AB)
					V=self.Vij(r, self.eps_AB, self.sig_AB, self.ron_AB, self.rcut_AB)
				temp=rvec*Vp
				energia+=V
				gradVector[i] +=  temp 
				gradVector[j] += -temp                
		return energia/Natoms,gradVector/Natoms


	# Calculate EXTENSIVE Energy and Gradient of the smoothened LJ potential
	def CalculateEnergySlower(self,snapshot):
		energia=0
		Natoms=snapshot.particles.N
		L=snapshot.box.Lx #Assume cubic box
		for i in range(Natoms):
			for j in range(i+1,Natoms):
				posi=np.array(snapshot.particles.position[i], dtype=np.float64)
				posj=np.array(snapshot.particles.position[j], dtype=np.float64)
				typei=snapshot.particles.typeid[i]
				typej=snapshot.particles.typeid[j]
				rvec=med.PeriodicDisplacement(posi,posj,L)
				r=np.sqrt(np.square(rvec).sum())

				#See what group the particles belong to
				if typei==typej:
					if typej==0:
						if r>self.rcut_AA:
							continue
						V=self.Vij(r, self.eps_AA, self.sig_AA, self.ron_AA, self.rcut_AA)
					else:
						if r>self.rcut_BB:
							continue
						V=self.Vij(r, self.eps_BB, self.sig_BB, self.ron_BB, self.rcut_BB)
				else:
					if r>self.rcut_AB:
						continue
					V=self.Vij(r, self.eps_AB, self.sig_AB, self.ron_AB, self.rcut_AB)
				energia+=V
		return energia




class JammingSphere:
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






#############################
# Old Potential definitions #
#############################

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

