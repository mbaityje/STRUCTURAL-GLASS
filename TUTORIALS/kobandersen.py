#!/usr/bin/python
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md


#Returns a hoomd.md.pair with Kob-Anderson Binary Lennard-Jones parameters
def KApotential(NeighborsListLJ):
    r_cutoff=2.5
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

