## module_phop.py     # this file's name #
## Iterative Barycenters Separation (IBS)  algorithm designed by Candelier & Dauchot.
## compute phop down to rather low phop values, i.e. use a low threshold, 
##      and then record all the phops of each event.
## 
## Francois P. Landes 2016 
##
from __future__ import print_function
import numpy as np
from numba import jit
from numba import float32, int32, int64, float64

import module_PBC

## always check PBCs as much as you can... don't foret them !!!

############################
#TODO: if memory is too big for our RAM, cut Natoms in independent pieces selctedAtoms
# and allow for the script to be called several time??? Or simply call it once but load pieces of trajectory one at a time 
# which would be more convenient for the writin toa gsd at the end.




@jit(nogil=True, nopython=True)
def sum_vector(vector):
    return vector[0]+vector[1]+vector[2]

@jit(nogil=True, nopython=True)
def xi_tc_T(tc,T):
    return abs(tc*1.0/(T-1.0) * (1.0-tc*1.0/(T-1.0)))**0.5
    # note that we use T-1 because in th eppaer it is between [0,T] and tehy use T: 
    #3 but in their case T itself is included, not in our C++-notations.
    ## thus our length is actually TB-TA-1 (number of steps studied is TB-TA-1+1)

@jit(nogil=True, nopython=True)
def unrolled_traj_func(trajectory,Lx, atom):
    trajDuration=len(trajectory)
    displacements = np.zeros((trajDuration, 3)) #, dtype=np.float)
    for t0 in range(1,trajDuration):
        displacements[t0] = module_PBC.applyPBC_vector_func(trajectory[t0,atom,:]-trajectory[t0-1,atom,:], Lx)    
    previousPosition = trajectory[0,atom,:]
    unrolled_traj = trajectory[:,atom,:].copy()
    for t0 in range(1,trajDuration):
        unrolled_traj[t0] = unrolled_traj[t0-1] + displacements[t0]
    return unrolled_traj
    
    
    
@jit(nogil=True, nopython=True)
def cumsum_func(vector, displacements, Lx): ## same as np.cumsum_func(vector) but in numba 
    N = len(vector)
    dimensionality=3
    cumsumA   = np.zeros((N,dimensionality))
    cumsumA_2 = np.zeros((N,dimensionality))
    
    ## TODO : optimization : doing the loop on each componentn 0,1,2 : is it faster ?
    cumsumA  [0,:] = vector[0,:]
    cumsumA_2[0,:] = vector[0,:]**2
    previousPosition = vector[0,:]
    for t0 in range(1,N):
#        displacement = displacements[t0] # ==  module_PBC.applyPBC_vector_func(vector[t0]-vector[t0-1], Lx)
        correctedPosition = previousPosition + displacements[t0]
        previousPosition = correctedPosition
        cumsumA  [t0,:] = cumsumA  [t0-1,:] + correctedPosition
        cumsumA_2[t0,:] = cumsumA_2[t0-1,:] + correctedPosition**2
    #print(previousPosition , vector[N-1])  #  -> is ofund to be 0 or +/- Lx 
    return (cumsumA, cumsumA_2)

#### unusued function #####
# @ jit(nogil=True, nopython=True)
#def reverse_vector(vector):
#    N = len(vector)
#    dimensionality=3
#    rev = np.zeros((N,dimensionality))
#    for t0 in range(0,N):
#        rev[t0,:] = vector[N-1-t0,:]
#    return rev

# not working ??? YES WORKING !! 
## TODO: check again and erase all refs saying it  does not work !! 
@jit(nogil=True, nopython=True)
def reverse_cumsum_func(cumsumX): ## same as np.cumsum(cumsumX[::-1]) but in numba
    N = len(cumsumX)
    dimensionality=3
    reversed_cumsumX = np.zeros((N,dimensionality))
    summ = cumsumX[N-1]
    for t0 in range(1,N):
        reversed_cumsumX[t0-1] = summ - cumsumX[N-1-t0]
    reversed_cumsumX[N-1] = summ 
    return reversed_cumsumX


## main function: returns phop[tc] with a range [TA,TB]
@jit(nogil=True , nopython=True)
def phopValue_func(subTraj, displacements, TA, TB, Lx):

    cumsumA, cumsumA_2 = cumsum_func(subTraj, displacements, Lx)
    cumsumB   = reverse_cumsum_func(cumsumA  )# not working  (not reversible)
    cumsumB_2 = reverse_cumsum_func(cumsumA_2)# not working
    ## slower method:
    ##rev = reverse_vector(subTraj)
    ##cumsumB, cumsumB_2 = cumsum_func(rev , Lx)

    phop = np.zeros(TB-TA)
    for tc in range(TA,TB,1):
        # indices LOCal to the arrays, for which time starts starts 0:
        tcLocA = tc-TA 
        tcLocB = TB-tc-1
        
        ## here we don't have to sum over t1 and t2 as in the Naive code version,
        ## which means that we do not compute many avergages, we just computed 
        ## the largest one once
        ## (+saved all intermediate results in the cumsum thingies)
        ## here we used that (a-b)^2 = a^2 - 2ab + b^2, with a=trajectory[t2,i]
        ## and b=<trajectory[tc:TB]> (see naive implementation for comparison)

        ## are these tmp variable really saving time?:
        averageA = cumsumA[tcLocA]/(tcLocA+1.0) 
        averageB = cumsumB[tcLocB]/(tcLocB+1.0)
        d1 = sum_vector( \
            cumsumA_2[tcLocA]/(tcLocA+1.0) \
            - 2.0 * averageA  * averageB \
            + (averageB)**2 )
        d2 = sum_vector( \
            cumsumB_2[tcLocB]/(tcLocB+1.0) \
            - 2.0 * averageB * averageA  \
            + (averageA )**2 )
        phop[tcLocA] = (d1*d2)**0.5 * xi_tc_T(tcLocA,TB-TA)
        
    tcLocA = np.argmax(phop)
    averageA = cumsumA[tcLocA]/(tcLocA+1.0)        # 3d vector
    Prev_StdDev = sum_vector(cumsumA_2[tcLocA]/(tcLocA+1.0) - averageA**2.0)   #scalar

    tcLocB = TB-TA-tcLocA-1    # TB-(tc-TA)-TA-1
    averageB = cumsumB[tcLocB]/(tcLocB+1.0)        # 3d vector
    Next_StdDev = sum_vector(cumsumB_2[tcLocB]/(tcLocB+1.0) - averageB**2.0)   #scalar
    
    average_distAB = (  sum_vector( (module_PBC.applyPBC_vector_func(averageA-averageB,Lx))**2 )  )**0.5

    ## dist. between 2 averages: alreadt computed in phop. 
    return phop, averageA, Prev_StdDev, Next_StdDev, average_distAB



## other algorithm, with sliding window (non optimal implementation, most likely, but whatever)
@jit(nogil=True , nopython=True)
def FixedWindow_phopValue_func(subTraj, displacements, TA, TB, Lx):
    ## tc-TA == TB-tc == Ttot/2 ##
    
    cumsumA, cumsumA_2 = cumsum_func(subTraj, displacements, Lx)
    cumsumB   = reverse_cumsum_func(cumsumA  )
    cumsumB_2 = reverse_cumsum_func(cumsumA_2)
    ## we compute twice as much as needed, but it's ok because we don't like this algorithm

    phop = np.zeros(TB-TA)
    tc = (TA+TB)/2
    # indices LOCal to the arrays, for which time starts starts 0:
    tcLocA = tc-TA 
    tcLocB = TB-tc-1

    ## are these tmp variable really saving time?:
    averageA = cumsumA[tcLocA]/(tcLocA+1.0) 
    averageB = cumsumB[tcLocB]/(tcLocB+1.0)
    d1 = sum_vector( \
        cumsumA_2[tcLocA]/(tcLocA+1.0) \
        - 2.0 * averageA  * averageB \
        + (averageB)**2 )
    d2 = sum_vector( \
        cumsumB_2[tcLocB]/(tcLocB+1.0) \
        - 2.0 * averageB * averageA  \
        + (averageA )**2 )
    phop[tcLocA] = (d1*d2)**0.5 # * xi_tc_T(tcLocA,TB-TA)  ## XXX TODO ! here I removed the xi, as they did in their papers !
    
    Prev_StdDev = sum_vector(cumsumA_2[tcLocA]/(tcLocA+1.0) - averageA**2.0)   #scalar
    Next_StdDev = sum_vector(cumsumB_2[tcLocB]/(tcLocB+1.0) - averageB**2.0)   #scalar
    average_distAB = (  sum_vector( (module_PBC.applyPBC_vector_func(averageA-averageB,Lx))**2 )  )**0.5

    ## dist. between 2 averages: alreadt computed in phop. 
    return phop, averageA, Prev_StdDev, Next_StdDev, average_distAB
################################################################################


################################################################################


## IBS ##
@jit(nogil=True, nopython=True)
def phopFromTrajectory_func(trajectory, Lx, p_threshold_low, every_recAlways):
    trajDuration = len(trajectory)
    Natoms = len(trajectory[0])
    ## mostly used for our own curiosity, not a final feature:
#    p_hop_hierarchy = [ [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] ]
    N_ObseravblesRecorded=11
    ## we set the value to the flag -42 so that it is ibvious a box is empty
    ## when the bix has not been filled with the correct value
    resultsArray = np.zeros((Natoms, trajDuration, N_ObseravblesRecorded), dtype=float32) #- 42.0
    Nact_atom    = np.zeros((Natoms), dtype=int32)
    hop_number=0
    OldIntervals_lengthMax = int(3e7)
    NewIntervals = np.zeros( (OldIntervals_lengthMax, 2) , dtype=int32) - 42
    OldIntervals = np.zeros( (OldIntervals_lengthMax, 2) , dtype=int32) - 42
    OldIntervals_lengthUsed_max = 0
    
    for atom in range(Natoms): # Natoms , but then we will do it in batches with numpy
        OldIntervals_lengthUsed = 0
        OldIntervals[OldIntervals_lengthUsed,0] = 0
        OldIntervals[OldIntervals_lengthUsed,1] = trajDuration 
        OldIntervals_lengthUsed = 1
#        OldIntervals = [[0, trajDuration]]
#        p_hop_hierarchy_depth=-1
        NactiveTimes=1

        displacements = np.zeros((trajDuration, 3)) #, dtype=np.float)
        for t0 in range(1,trajDuration):
            displacements[t0] = module_PBC.applyPBC_vector_func(trajectory[t0,atom,:]-trajectory[t0-1,atom,:], Lx)

        while OldIntervals_lengthUsed > 0: #len(OldIntervals) > 0 :
#            p_hop_hierarchy_depth+=1
            NewIntervals_lengthUsed=0  #NewIntervals=[]
            # TA and TB are indeed indices, i.e. starting from 0 and ending at Nframes-1.
#            print(OldIntervals)
            for TATBindex in range(OldIntervals_lengthUsed): #for (TA,TB) in OldIntervals:
                TA = OldIntervals[TATBindex,0]
                TB = OldIntervals[TATBindex,1]
                ## dismiss the search when two consecutive "events" have been found: 
                if TB-TA>1: 
                    subTraj = trajectory[TA:TB,atom]

                    ##### core fo the algorithm: ###########
                    phop, averageA, Prev_StdDev, Next_StdDev, average_distAB = phopValue_func(subTraj, displacements[TA:TB], TA, TB, Lx)
                    ########################################
                    tcLocA = np.argmax(phop)
                    p_hop_value = phop[tcLocA]
                    tc = tcLocA+TA
                    
#                    if p_hop_value > StdDevInternalA+StdDevInternalB:        
                    
                    if p_hop_value > p_threshold_low:
#                        p_hop_hierarchy[p_hop_hierarchy_depth].append(p_hop_value)
                        resultsArray[atom, tc, 0] = tc
                        resultsArray[atom, tc, 1] = p_hop_value
                        resultsArray[atom, tc, 2] = 1
                        resultsArray[atom, tc, 3] = TA
                        resultsArray[atom, tc, 4] = TB
                        resultsArray[atom, tc, 5] = averageA[0]
                        resultsArray[atom, tc, 6] = averageA[1]
                        resultsArray[atom, tc, 7] = averageA[2]
                        resultsArray[atom, tc, 8] = Prev_StdDev
                        resultsArray[atom, tc, 9] = Next_StdDev
                        resultsArray[atom, tc,10] = average_distAB
#                        resultsArray[atom, NactiveTimes, 11] = tc
                        NactiveTimes+=1  ## just for that atom in particular !!
                        
                        ## TODO : after all events for all atoms have been recorded to this, 
                        ## optimize file length (compress) bu recording inly the list of events' properties for eaxh atom.
                        ## the output file would then be shape  (Natoms, number-of-events-of-most-active-atom)
                        ## + the positions of full system for each active time (ooking at unique times)
                        ## XXX however this may turn out to be useless if I'm lucky.
                    
                        hop_number+=1

                        ## uodate of the "list" of intervales to search events in ##                        
                        ###NewIntervals.append([TA,tc])
                        ###NewIntervals.append([tc,TB])
                        NewIntervals[NewIntervals_lengthUsed  ,0] = TA
                        NewIntervals[NewIntervals_lengthUsed  ,1] = tc 
                        NewIntervals[NewIntervals_lengthUsed+1,0] = tc
                        NewIntervals[NewIntervals_lengthUsed+1,1] = TB 
                        NewIntervals_lengthUsed += 2
                        if NewIntervals_lengthUsed >= OldIntervals_lengthMax:
                            print('error: increase size of NewIntervals array to avoid this error')
                            return resultsArray, Nact_atom, OldIntervals_lengthUsed_max, hop_number
                    else: 
                        ## there is no hop detected in that interval [TA,TB]:
                        for tc in range(TA,TB):
                            if tc%every_recAlways==0:
                                ## we want to record frames anyways, every some...
                                resultsArray[atom, tc, 0] = tc
                                resultsArray[atom, tc, 1] = -1.0*phop[tc-TA] 
                                resultsArray[atom, tc, 2] = 0
                                resultsArray[atom, tc, 3] = TA
                                resultsArray[atom, tc, 4] = TB
                                resultsArray[atom, tc, 5] = 42
                                resultsArray[atom, tc, 6] = 42
                                resultsArray[atom, tc, 7] = 42
                                resultsArray[atom, tc, 8] = -1
                                resultsArray[atom, tc, 9] = -1
                                resultsArray[atom, tc,10] = -1
                                ## resultsArray[atom, tc, 11] is not used because the site is not in the lsit of active times ##
                                
            #### end of the for loop on all (TA,TB) pairs ###
            OldIntervals_lengthUsed = NewIntervals_lengthUsed   
            OldIntervals_lengthUsed_max = max(OldIntervals_lengthUsed_max, OldIntervals_lengthUsed) 
            OldIntervals[:OldIntervals_lengthUsed] = NewIntervals[:NewIntervals_lengthUsed] ## almost a ffull copy !
        
#        ## remember the size of the "list"
#        resultsArray[atom, 1:NactiveTimes,11] = np.sort(resultsArray[atom,1:NactiveTimes,11])
#        resultsArray[atom, 0 ,11] = NactiveTimes
        Nact_atom[atom] = NactiveTimes-1
        
    ###############################################################################
    return resultsArray, Nact_atom, OldIntervals_lengthUsed_max, hop_number


################################################################################

## sliding window algorithm ##
@jit(nogil=True, nopython=True)
def FixedWindow_phopFromTrajectory_func(trajectory, Lx, p_threshold_low, every_recAlways):
    trajDuration = len(trajectory)
    Natoms = len(trajectory[0])
#    uniqueHops=[]
    ## mostly used for our own curiosity, not a final feature:
#    p_hop_hierarchy = [ [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] ]
    N_ObseravblesRecorded=12
    ## we set the value to the flag -42 so that it is ibvious a box is empty
    ## when the bix has not been filled with the correct value
    resultsArray=np.zeros((Natoms, trajDuration, N_ObseravblesRecorded), dtype=float32) #- 42.0
    Nact_atom    = np.zeros((Natoms), dtype=int32)
    hop_number=0
    
    ## TODO: make this width an argument of this function, maybe ?
    TF=5
    
    
    for atom in range(Natoms): # Natoms , but then we will do it in batches with numpy
        NactiveTimes=1
        displacements = np.zeros((trajDuration, 3)) #, dtype=np.float)
        for t0 in range(1,trajDuration):
            displacements[t0] = module_PBC.applyPBC_vector_func(trajectory[t0,atom,:]-trajectory[t0-1,atom,:], Lx)
            TA = max(t0-TF,0)
            TB = min(t0+TF,trajDuration-1)
            subTraj = trajectory[TA:TB,atom]

            ##### core fo the algorithm: ###########
            phop, averageA, Prev_StdDev, Next_StdDev, average_distAB = FixedWindow_phopValue_func(subTraj, displacements[TA:TB], TA, TB, Lx)
            ########################################
            tcLocA = t0-TA #np.argmax(phop)
            p_hop_value = phop[tcLocA]
            tc = t0 #  tcLocA+TA
            
            if p_hop_value > p_threshold_low:
                resultsArray[atom, tc, 0] = tc
                resultsArray[atom, tc, 1] = p_hop_value
                resultsArray[atom, tc, 2] = 1
                resultsArray[atom, tc, 3] = TA
                resultsArray[atom, tc, 4] = TB
                resultsArray[atom, tc, 5] = averageA[0]
                resultsArray[atom, tc, 6] = averageA[1]
                resultsArray[atom, tc, 7] = averageA[2]
                resultsArray[atom, tc, 8] = Prev_StdDev
                resultsArray[atom, tc, 9] = Next_StdDev
                resultsArray[atom, tc,10] = average_distAB
#                resultsArray[atom, NactiveTimes, 11] = tc
                NactiveTimes+=1  ## just for that atom in particular !!
                hop_number+=1
                
#        ## remember the size of the "list"
#        resultsArray[atom, 1:NactiveTimes,11] = np.sort(resultsArray[atom,1:NactiveTimes,11])
#        resultsArray[atom, 0 ,11] = NactiveTimes
        Nact_atom[atom] = NactiveTimes-1
    ###############################################################################
    return resultsArray, Nact_atom, hop_number








#maxNumberOfHops = trajDuration
## useless if maxNumberOfHops == trajDuration  :
### if you need to use less memory and expect few events:
##maxNumberOfHops = 1000 
#maxNumberOfHops_growthFactor = 2 


### stupid function, mostly useless, replaces "np.resize"  which seems broken 
### for multi -dimensional arrrays
#def extend_func(resultsArray, maxNumberOfHops_growthFactor):
#    maxNumberOfHops = np.shape(resultsArray)[1]
#    new_maxNumberOfHops = int( maxNumberOfHops * maxNumberOfHops_growthFactor )
#    new_resultsArray = np.zeros((Natoms, new_maxNumberOfHops, N_ObseravblesRecorded))-42
#    new_resultsArray[:,:maxNumberOfHops,:] = resultsArray
#    maxNumberOfHops = new_maxNumberOfHops
#    print 'resizing resul array to accomodate more phops: it now goes up to: ', maxNumberOfHops, ' phops (for some atoms only of course)'
#    resultsArray = new_resultsArray
#    del new_resultsArray
#    return resultsArray, maxNumberOfHops
###if hop_number >= maxNumberOfHops :
###    resultsArray, maxNumberOfHops = extend_func(resultsArray, maxNumberOfHops_growthFactor)

################################################################################









##    ### DEBUG: (INCORRECT) Naive implementation (wothout the PBC !!!)
#        d1good=d1
#        d2good=d2
##    phopNaive= np.zeros(TB-TA)
##    for tc in range(TA,TB,1):
#        d1=0.0
#        for t1 in range(TA,tc+1): 
#            d1 += sum_vector( (trajectory[t1,i] - np.mean(trajectory[tc:TB,i] ,0) )**2 )
#        d1 /= (tc-TA+1)
#        d2=0.0
#        for t2 in range(tc, TB):
#            d2 += sum_vector( (trajectory[t2,i] - np.mean(trajectory[TA:tc+1,i] ,0) )**2 )
#        d2 /= (TB-tc)
#        phopNaive[tc-TA ] = (d1*d2)**0.5  *xi[tcLocA]
#    print( phop - phopNaive)
#[ -1.52655666e-16  -3.27862737e-16  -4.18068358e-16  -3.92047506e-16  -1.24900090e-15  -4.33680869e-16]
# end of DEBUG: end of the naive implementation, to do a gloabl test ##


## DEBUG
#a-phop[::-1]
#Out[220]: 
#array([  0.00000000e+00,  -8.67361738e-19,   0.00000000e+00,
#         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
#         0.00000000e+00,   0.00000000e+00,   8.67361738e-19,
#         0.00000000e+00])

################################################################################

#print 'finished searching for events. '
### some points may have been missed by the algorithm (I do
### not fully undestand why, but..)
#for atom in [1]: #range(Natoms):
#    for tc in range(0,trajDuration, every_recAlways):
#        if resultsArray[atom,tc,0]< -1 :
#            TA=tc
#            while resultsArray[atom,TA,0]< -1 and TA>0:
#                TA-=1
#            TB=tc
#            while resultsArray[atom,TB,0]< -1 and TB<trajDuration-1:
#                TB+=1
#            print 'Weird point found in atom ', atom, ' index tc=', tc
#            if TB-TA>1:
#                subTraj = trajectory[TA:TB, atom]
#                phop = phopValue_func(subTraj, TA, TB)
#                
#                for tc_around in range(TA,TB,1):
#                    if tc_around%every_recAlways==0:
#                        resultsArray[atom, tc_around, 0] = tc
#                        resultsArray[atom, tc_around, 1] = phop[tc_around-TA] # tcLocA = tc-TA
#                        resultsArray[atom, tc_around, 2] = 0
#                        resultsArray[atom, tc_around, 3] = TA
#                        resultsArray[atom, tc_around, 4] = TB
#                
#            
    
### malfunctioning piece of code: replaced with simpler and maybe slightly less efficient code, using simply tc%every_recAlways==0 istead of THAT:
#TAprime = ( TA/every_recAlways + 1) * every_recAlways
#if TA%every_recAlways==0:
#    TAprime -= every_recAlways
#TBprime = TB-(TB%every_recAlways)   
#if TAprime<0:
#    TAprime=0
#if TBprime>trajDuration-1:
#    TBprime=trajDuration-1
#explore = range(TAprime, TBprime+1, every_recAlways) #+[TAprime]
#for tc_around in explore:
#    resultsArray[atom, tc_around, 0] = tc
#    resultsArray[atom, tc_around, 1] = phop[tc_around-TAprime] # tcLocA = tc-TA
#    resultsArray[atom, tc_around, 2] = 0
#    resultsArray[atom, tc_around, 3] = TA
#    resultsArray[atom, tc_around, 4] = TB

################################################################################

### add a frame every  "every_recAlways"  so that you can also have 
### regualrly-spaced records of what happened
#tcs = np.insert(tcs,np.arange(0,trajDuration,every_recAlways),0)

