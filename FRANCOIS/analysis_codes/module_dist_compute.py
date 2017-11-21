#module_dist_compute.py
from __future__ import print_function
import numpy as np
from numba import float32, int32, int64, float64
from numba import jit
import numba
import module_PBC
# tricks about numba:
# 1) comment out all the @jit to debug 
# 2) if @jit(options) does not work, try  @jit() alone
# 3) prints statements are a mess

## this is a container for the core @jit function
## it does not need speed.
## TODO : make the beggining of this (numba)jit-compatible, so as to be able to call this function from a numba function
def dist_computer(positions, Lx, expectedDensity, shift_vector, rMax=-1):
    '''
    computes all the pair distances for all the atoms present in "positions". 
    It also provides the atomIds of each pair.
    It limits the computation to r<rMax.
    It is fast thanks to cell decomposition (but no parallel computing)
    '''
    verbose=0
    
    if rMax == -1: 
        rMax=3.0
        rMax=Lx/2.0 
        ## Lx/2 seems like the most reasonable value, becasue is the largest allowed .. right ?!

    if rMax > Lx/2.0:
        rMax = Lx/2.0
        #rMax = np.round(rMax,1)
        print(" if you take  rMax  more than  Lx/2.0  you are going to count some neighbors as if they were far from you, when actually they are close to you !\nWe correct rMax to its maximal value of Lx/2=", rMax)

    ####################################################
    ######## unphyscial, irrelevant choices: ###########

    ## actual number of cells is   linear_number_cells^3 
    ## heuristically, the optimal speedup ??
    linear_number_cells = int(2.0 * Lx / (rMax/2.0)) 
    size_factor = 1
    shift_cells = 0
    if shift_vector[0] > 0: 
        size_factor = int(Lx / shift_vector[0] + 1e-9)
        linear_number_cells -= (linear_number_cells % size_factor)    ## because we then multiply by this factor factor_size
        if linear_number_cells <= 1:
            linear_number_cells = size_factor
        shift_cells = linear_number_cells / size_factor     ## this MUST be an exact division
        if shift_cells % 1 > 1e-9:
            print("\nwrong!! The linear_number_cells should be an exact multiple of size_factor !!\n\n")
    if verbose:
        print('using linear_number_cells=', linear_number_cells, ' , i.e. a total number of cells of ', linear_number_cells**3)
    # then we will do:  cellsUnitLength = Lx*1.0/linear_number_cells

#    print("shift_cells, shift_vector, ...", shift_cells, shift_vector)

    #####################################################
    ######### beginning of the "main" function  #########

    #### DEBUG: you may want to do for smaller systems sometimes
    #    Natoms=6400
    #    positions=positions[:Natoms]
    Natoms = positions.shape[0]
    if verbose:
        print('Natoms='+str(Natoms))

    ## length of each cells in units of sigma_AA=1
    cellsUnitLength = Lx*1.0/linear_number_cells

    ###################################################
    ##### prepare some arrays, using pure python: #####

    ## define the mesh or cell list into which we sort the particles
    cellsGrid = np.arange(-Lx/2,Lx/2+cellsUnitLength, cellsUnitLength)
    if np.min(positions)< np.min(cellsGrid) or np.max(positions) > np.max(cellsGrid) :
        print('error: your cellsGrid does not contain all particles. \nExit.\n')
        raise SystemExit

    ## this code currently supports only 3, because there are 3 loops 
    ## inside the bulk of code: i,j,k (self) and then i,j,k (neighbors)
    dim = 3 
    max_numberDensity = 4.0*expectedDensity # unphysical, need to be larger than actual number density.
    max_occupancy     = int(max_numberDensity*cellsUnitLength**dim)
    ### XXX memory errors may come from this: (it needs to be huge because we have so many enighbors
    ###  TODO: take the 2nd colum of this resultArray seoparately, declare it with type uint16 or uint32 to save some memory
    max_neighb_number = int(max_numberDensity*5*rMax**dim)
    
    back_table = np.zeros(((linear_number_cells,)*dim+(max_occupancy,)), dtype=np.int)
    nOcc_array = np.zeros((linear_number_cells,)*dim, dtype=np.int) # keep track of how full the local column of back_table[i,j,k,:] is.



    #########################################################
    ##### the core, the numba function: #####################
    #@numba.jit( float64[:](float64[:,:], int64[:,:,:], int64[:,:,:,:],   float64, float64, int64, float64))
    @jit(nopython=True, nogil=True)
    def dist_numbaFunction(positions, nOcc_array, back_table, rMax, linear_number_cells, cellsUnitLength, shift_vector):
    
        Natoms = positions.shape[0]
        resultArray = np.zeros((Natoms,max_neighb_number,2))
        dist_num = np.zeros(Natoms, dtype=int32)
        
        ## sweep atom once to build their index:
        histog  = np.digitize(positions, cellsGrid) 
        histog -= 1     ## becasue I don't like numpy's convention
        ### sweep atoms once to build the backwards-table of index:
        for atom in range(Natoms):
            i = histog[atom,0]
            j = histog[atom,1]
            k = histog[atom,2]
            back_table[i,j,k,nOcc_array[i,j,k]] = atom
            nOcc_array[i,j,k] += 1
#            ## We want to check that, here, but numba does not agree:
#            if nOcc_array[i,j,k] == max_occupancy :
#                print("Error: increase max_numberDensity so that max_occupancy is larger than "+str(max_occupancy))
#                raise SystemExit

        limitSup = (linear_number_cells+3**0.5) # # Lx/cellsUnitLength + 3**0.5 +1

        indexMax = int(np.ceil( (rMax/ cellsUnitLength) +1 ))
#        print('indexMax=')
#        print(indexMax)
        # max index shift (in units of cell size) where we will look for neighbors (at distance less than rMax)

        NMax = 8000
        linear_number_cellsStart = linear_number_cells
        if linear_number_cellsStart**3*1.2 > NMax :
            linear_number_cellsStart = int((NMax/1.2)**0.33333333333333333333333)
            

        previousRec = 0.0
        ### sweep each atom once and then all of its neighbors (some of them are uselss, not many)
        for i in range(linear_number_cellsStart):
            for j in range(linear_number_cellsStart):
                for k in range(linear_number_cellsStart):
                    ## now looking for neighbors of this cell
                  
                    ## TODO : optimization : 
                    ## loop only over shorter ranges, that may go outsuide the box,
                    ## and then use  ineighb % linear_number_cells  to go back inside the box
                    # the ranges may be sthg like : range(i, i+indexMax):   #or more cost: (i-indexMax, i+indexMax) (then indexMax must not be too large)
                    for ineighb in range(linear_number_cells): 
                        for jneighb in range(linear_number_cells):
                            for kneighb in range(linear_number_cells):
                                
                                ## super-optimization: we cut the search to the enclosing sphere.
                                ## apply PBC to the *distances* 
                                di = module_PBC.applyPBC_pair_func(i, ineighb, linear_number_cells)
                                if shift_vector[0]>0:
                                    di = module_PBC.applyPBC_pair_func(i, ineighb - shift_cells, linear_number_cells) # XXX
                                
                                
                                dj = module_PBC.applyPBC_pair_func(j, jneighb, linear_number_cells)
                                dk = module_PBC.applyPBC_pair_func(k, kneighb, linear_number_cells)
#                                di = abs(i-ineighb) % linear_number_cells/2    ## rounding down is the desired behav
#                                dj = abs(j-jneighb) % linear_number_cells/2 
#                                dk = abs(k-kneighb) % linear_number_cells/2
                                if np.sqrt((di)**2+(dj)**2+(dk)**2) <= limitSup: 

                                    ### this could be placed earlier, but here is fine too (faster):
                                    nOccLoc = nOcc_array[i,j,k]
                                    nOccNeighb = nOcc_array[ineighb,jneighb,kneighb]
                                    
                                    for nSelf in range(nOccLoc):
                                        atom = back_table[i,j,k,nSelf]
                                        pos = positions[atom]
                                        for nNeighb in range(nOccNeighb):
                                            atom_b = back_table[ineighb,jneighb,kneighb,nNeighb]
                                            pos_neighb = positions[atom_b]
                                            displacements = module_PBC.applyPBC_vector_func(pos - pos_neighb - shift_vector, Lx) 
                                            dist = np.sqrt( np.sum((displacements)**2) )    ## TODO optimize this as a separate function,
                                                        # knowing that dim=3
                                            
                                            if dist < rMax :
                                                resultArray[atom, dist_num[atom],0] = dist
                                                resultArray[atom, dist_num[atom],1] = atom_b
                                                dist_num[atom] += 1
                                                
                                                ##DEBUG : used to find and check the limitSup :
                                                #tmp= np.sqrt((di)**2+(dj)**2+(dk)**2)
                                                #if tmp > previousRec :
                                                #   previousRec = tmp
                                                #   save =((di), (dj), (dk))
                                                
        #### DEBUG too:
        #print(previousRec, limitSup , rMax/cellsUnitLength, cellsUnitLength, rMax)
        #print(save)
        return resultArray, dist_num
        
    ############################################################################
    ################ CORE FUNCTION (CALL): #####################################
    resultArray, dist_num = dist_numbaFunction(positions, nOcc_array, back_table,  rMax, linear_number_cells, cellsUnitLength, shift_vector)
    ############################################################################
    
    return resultArray, dist_num
                             
