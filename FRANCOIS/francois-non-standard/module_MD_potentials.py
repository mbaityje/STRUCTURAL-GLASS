## module containg the definition and parameters and smoothing parameters of the potentials you may want to use.
## list:
# LJ
# WCA 

from hoomd import md

def set_interaction_potential_LJ(NeighborsListLJ):
    r_cutoff=2.5
    eps_AA=1 
    eps_AB=1.5
    eps_BB=0.5
    sig_AA=1
    sig_AB=0.8
    sig_BB=0.88
    r_on_cutoff=1.2
    ## specify Lennard-Jones interactions between particle pairs
    myLjPair = md.pair.lj(r_cut=r_cutoff, nlist=NeighborsListLJ)
    myLjPair.pair_coeff.set('A', 'A', epsilon=eps_AA, sigma=sig_AA, r_cut=r_cutoff*sig_AA, r_on=r_on_cutoff*sig_AA)
    myLjPair.pair_coeff.set('A', 'B', epsilon=eps_AB, sigma=sig_AB, r_cut=r_cutoff*sig_AB, r_on=r_on_cutoff*sig_AB)
    myLjPair.pair_coeff.set('B', 'B', epsilon=eps_BB, sigma=sig_BB, r_cut=r_cutoff*sig_BB, r_on=r_on_cutoff*sig_BB)
    myLjPair.set_params(mode="xplor")   ## smooth interpolation starting at r_on and ending at r_cut
    return myLjPair

def WCA_smooth(r, rmin, rmax, epsilon, sigma):
    V = 4 * epsilon * ( (sigma / r)**12 - (sigma / r)**6) + epsilon     -  epsilon * (36*2**(-1/3.)*(r-rmax)**2)     / sigma**2
    F = 4 * epsilon / r * ( 12 * (sigma / r)**12 - 6 * (sigma / r)**6)  +  epsilon * (36*2**(-1/3.)*(r-rmax)   ) *2  / sigma**2 
    return (V, F)

def set_interaction_potential_WCA(NeighborsListLJ):
    r_cutoff=2**(1./6)  ## 2**(1./6) == 1.122462048309373
    r_min = 0.722462048309373
    TABLEWIDTH=1000 ## number of particle pairs 
    
    eps_AA=1 
    eps_AB=1.5
    eps_BB=0.5
    sig_AA=1
    sig_AB=0.8
    sig_BB=0.88
    ## specify WCA interactions between particle pairs
    myLjPair = md.pair.table(width=TABLEWIDTH, nlist=NeighborsListLJ)
    myLjPair.pair_coeff.set('A', 'A', func=WCA_smooth, rmin=r_min*sig_AA, rmax=r_cutoff*sig_AA, coeff=dict(epsilon=eps_AA, sigma=sig_AA))
    myLjPair.pair_coeff.set('A', 'B', func=WCA_smooth, rmin=r_min*sig_AB, rmax=r_cutoff*sig_AB, coeff=dict(epsilon=eps_AB, sigma=sig_AB))
    myLjPair.pair_coeff.set('B', 'B', func=WCA_smooth, rmin=r_min*sig_BB, rmax=r_cutoff*sig_BB, coeff=dict(epsilon=eps_BB, sigma=sig_BB))
    return myLjPair

#def set_interaction_potential_WCA(NeighborsListLJ):
#    r_cutoff=2**(1./6)
#    eps_AA=1 
#    eps_AB=1.5
#    eps_BB=0.5
#    sig_AA=1
#    sig_AB=0.8
#    sig_BB=0.88
#    ## specify WCA interactions between particle pairs
#    myLjPair = md.pair.lj(r_cut=r_cutoff, nlist=NeighborsListLJ)
#    myLjPair.pair_coeff.set('A', 'A', epsilon=eps_AA, sigma=sig_AA, r_cut=r_cutoff*sig_AA)
#    myLjPair.pair_coeff.set('A', 'B', epsilon=eps_AB, sigma=sig_AB, r_cut=r_cutoff*sig_AB)
#    myLjPair.pair_coeff.set('B', 'B', epsilon=eps_BB, sigma=sig_BB, r_cut=r_cutoff*sig_BB)
#    myLjPair.set_params(mode="shift")
#    return myLjPair
