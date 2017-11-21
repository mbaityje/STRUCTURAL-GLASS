## module_descriptors.py : 
## called from:
# - construct_TrainingSet.py  
# - assign_softness.py 
## It contains the core function(s) for computing descriptors, given the neighboring positions.
## it makes use of numba for performance reasons.
from __future__ import print_function

import numpy as np
import numba 
from numba import float32, int32, int64, float64
from numba import jit

import module_PBC
import module_dist_compute


## TODO: code for other families of descriptor functions here. 
@jit(nopython=True, nogil=True)
def descriptor_func_singleDistance(descriptor_family_name, descriptor_family_parameters, dist):
    if descriptor_family_name == 1:
        return np.exp( -(dist-descriptor_family_parameters[:,0])**2/descriptor_family_parameters[:,1]**2)
    if descriptor_family_name == 2 :
        pass 
        ## something else 

 
 ## function used only in: B1-construct_TrainingSet
@jit(nopython=True, nogil=True)
def descriptor_func_allPositions_oneSnapshot(resultArray, dist_num, descriptor_family_name, descriptor_family_parameters, atomTypes):
    
    Natoms = atomTypes.shape[0]
    Ntypes = int(np.max(atomTypes)+1)
    descriptors_size_m = len(descriptor_family_parameters[:,0]) # + the angular stuff
    # temporary arrays: #
    ## (we want extra precision because we will add many elements) ##
    all_atoms_descriptors_array = np.zeros((Natoms, descriptors_size_m*Ntypes))#, dtype=float32)
    descriptors_array           = np.zeros( (Ntypes, descriptors_size_m))#, dtype=float64) 
    descriptor_cutoff = np.max(descriptor_family_parameters[:,0]) - 1e-6


    for atom in range(Natoms):
        for num in range(dist_num[atom]):
            dist   = resultArray[atom, num ,0]
            atom_b = int(resultArray[atom, num ,1])
            # distance between 'atom' and 'atom_b' : dist 
            if dist < descriptor_cutoff :
                tmp = descriptor_func_singleDistance(descriptor_family_name, descriptor_family_parameters, dist)
            for i in range(descriptors_size_m):
                ## numba refuses to do my full copy :( using [:] .. it's stil faster !! :)
                descriptors_array[ atomTypes[atom_b],i ] += tmp[i]

        ## record the set of {G_i}'s for that instant
        for AtomType in range(Ntypes):
            all_atoms_descriptors_array[atom, descriptors_size_m*AtomType:descriptors_size_m*(AtomType+1)] = descriptors_array[AtomType]
        ## we do not need as much precision there (in the final result) ##
         
    return all_atoms_descriptors_array


 
 
 
## function used only in: B1-construct_TrainingSet
@jit(nopython=True, nogil=True)
def descriptor_func_allInstants(AtomTimesAll, trajectory, backArray, descriptor_family_name, descriptor_family_parameters, Lx, atomTypes):
    
    Natoms = len(trajectory[0,:,0])
#    Nexamples = len(AtomTimesAll)/2
    descriptors_size_m = len(descriptor_family_parameters[:,0])
    Ntypes = int(np.max(atomTypes)+1)
    all_instants_descriptors_array = np.zeros((len(AtomTimesAll), descriptors_size_m*Ntypes)) 

    ## choice of the cutoff: it will be assumed that the max of the 0th column of parmaeters gives a meanigful radius.    
    descriptor_cutoff = np.max(descriptor_family_parameters[:,0]) - 1e-6
    
    for instant_index in range(len(AtomTimesAll)): 
#        if instant_index % int(len(AtomTimesAll)/10) ==0:
#            print("1/10th done: instant is: ", int(instant_index))

        # temporary array: #
        ## (we want extra precision because we will add many elements) ##
        descriptors_array     = np.zeros( (Ntypes, descriptors_size_m), dtype=float64) 

        atom     = AtomTimesAll[instant_index,0]
        tc_index = backArray[AtomTimesAll[instant_index,1]]
        allpos   = trajectory[tc_index]
        pos      = allpos[atom]

        for atom_b in range(Natoms):
            # distance between 'atom' and 'atom_b' : dist 
            pos_neighb = allpos[atom_b]
            displacements = module_PBC.applyPBC_vector_func(pos - pos_neighb, Lx) 
            dist = np.sqrt( np.sum((displacements)**2) )    
            if dist < descriptor_cutoff :

                tmp = descriptor_func_singleDistance(descriptor_family_name, descriptor_family_parameters, dist)
                for i in range(descriptors_size_m):
                    ## numba refuses to do my full copy :( using [:] .. it's stil faster !! :)
                    descriptors_array[ atomTypes[atom_b],i ] += tmp[i]
        
        
        ## record the set of {G_i}'s for that instant
        for AtomType in range(Ntypes):
            all_instants_descriptors_array[instant_index, descriptors_size_m*AtomType:descriptors_size_m*(AtomType+1)] = descriptors_array[AtomType]
        ## we do not need as much precision there. ##
        
    return all_instants_descriptors_array

