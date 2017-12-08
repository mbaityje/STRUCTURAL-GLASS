#!/usr/bin/python
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md
import numpy as np
import module_measurements as med

class LJ():
    """ A class that contains all the information on the potential
    type: KA, KAshort

    """

    mode="xplor" #I only implement xplor, cause it's the only one I use. To obtain the other modes you only need to remove pieces and look at tutorials T11, T12.

    def __init__(self, NeighborsListLJ, type="KAshort"):
        if type=="KA": #Kob-Andersen potential
            self.setKA(NeighborsListLJ)
        elif type=="KAshort": #Kob-Andersen for small systems, as proposed by Heuer
            self.setKA(NeighborsListLJ, r_cutoff=1.4, r_buff=0.0)
        else:
            print("type Non implementato")
            raise SystemExit
 
    def setKA(self,NeighborsListLJ, r_on_cutoff=1.2, r_cutoff=2.5, r_buff=None):
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
            NeighborsListLJ.set_params(r_buff=self.r_buff)
        self.myLjPair = md.pair.lj(r_cut=self.r_cutoff, nlist=NeighborsListLJ)
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


    # Calculate Energy and Gradient of the smoothened LJ potential
    def CalculateGradient(self,snapshot):
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




#############################
# Old Potential definitions #
#############################





#Returns a hoomd.md.pair with Kob-Anderson Binary Lennard-Jones parameters
def KApotential(NeighborsListLJ, r_cutoff=2.5): 
    eps_AA=1 
    eps_AB=1.5
    eps_BB=0.5
    sig_AA=1
    sig_AB=0.8
    sig_BB=0.88
    r_on_cutoff=1.2
    myLjPair = md.pair.lj(r_cut=r_cutoff, nlist=NeighborsListLJ)
    myLjPair.pair_coeff.set('A', 'A', epsilon=eps_AA, sigma=sig_AA, r_cut=r_cutoff*sig_AA, r_on=r_on_cutoff*sig_AA)
    myLjPair.pair_coeff.set('A', 'B', epsilon=eps_AB, sigma=sig_AB, r_cut=r_cutoff*sig_AB, r_on=r_on_cutoff*sig_AB)
    myLjPair.pair_coeff.set('B', 'B', epsilon=eps_BB, sigma=sig_BB, r_cut=r_cutoff*sig_BB, r_on=r_on_cutoff*sig_BB)
    myLjPair.set_params(mode="xplor")
    return myLjPair

#Returns a hoomd.md.pair with Kob-Anderson Binary Lennard-Jones parameters
def KApotentialShort(NeighborsListLJ):
    eps_AA=1
    eps_AB=1.5
    eps_BB=0.5
    sig_AA=1
    sig_AB=0.8
    sig_BB=0.88
    r_on_cutoff=1.2
    r_cutoff=1.4
    NeighborsListLJ.set_params(r_buff=0.0)
    myLjPair = md.pair.lj(r_cut=r_cutoff, nlist=NeighborsListLJ)
    myLjPair.pair_coeff.set('A', 'A', epsilon=eps_AA, sigma=sig_AA, r_cut=r_cutoff*sig_AA, r_on=r_on_cutoff*sig_AA)
    myLjPair.pair_coeff.set('A', 'B', epsilon=eps_AB, sigma=sig_AB, r_cut=r_cutoff*sig_AB, r_on=r_on_cutoff*sig_AB)
    myLjPair.pair_coeff.set('B', 'B', epsilon=eps_BB, sigma=sig_BB, r_cut=r_cutoff*sig_BB, r_on=r_on_cutoff*sig_BB)
    myLjPair.set_params(mode="xplor")
    return myLjPair

