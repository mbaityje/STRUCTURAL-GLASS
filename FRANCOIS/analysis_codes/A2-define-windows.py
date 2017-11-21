## Here we just look at a light file describing the time of rearrangements (or no-rearrangements) and define the hard and soft (trainer) times accordingly. This is extremely light, as we do not look at the actual positions of any atom.
## -1- Get a given, large number of points, by autmatically computing the adequate slowWindow .... in order to study the sample' size dependence.

from __future__ import print_function
import gsd.pygsd
import gsd.hoomd
import gsd.fl

import sys
import numpy as np

import module_filenamesHandler 

######
## plots: 
import module_importPlotParams
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True          
plt.rcParams.update(module_importPlotParams.fparams())
plt.ion()   ## optional
######


print('usage:')
print('run  ',sys.argv[0],'   [filename]   [p_threshold_low=-1(autodetect)]   [p_threshold_high=(-)0.2]    [Nsamples=(-)100, slowWindow=?]   [fastWindow=2] ')
filename = sys.argv[1]
rootname = module_filenamesHandler.get_rootname(filename)
suffix   = module_filenamesHandler.get_suffix(filename)

p_threshold_low = -1
if len(sys.argv)>2:
    p_threshold_low = float(sys.argv[2])

p_threshold_high = -0.2
if len(sys.argv)>3:
    p_threshold_high = float(sys.argv[3])

slowWindow=-100
if len(sys.argv)>4:
    slowWindow = int(sys.argv[4])
if slowWindow < 0:
    Nsamples = - slowWindow
else:
    Nsamples = 42 # unused

fastWindow=2
if len(sys.argv)>5:
    fastWindow = int(sys.argv[5])

kind = 'B'
kind = 'A'
kind = 'all'
if len(sys.argv)>6:
    kind = sys.argv[6]



throwBoundaries_soft = 1

plotting=1

debug_or_choose_parameters = True
debug_or_choose_parameters = False


# above: from phop.py,may be useless.
## le read de decorated: issu de la partie plots de phop 
## ensuite: issu de la partie eventsData de phop, analyse qui a sa place ici.
## ensuite:partie defintion med soft hard de construct, qui doit etre ici et non la bas. 


################################################################################################################################################################
if "deco" in filename :   # read the decorated file produced before # 

    if p_threshold_low == -1 : 
        p_threshold_low = float(module_filenamesHandler.filename_parser(filename[:-4], 'pLow'))
        print('reading that the file has a p_threshold_low=', p_threshold_low, ' ...')
        if p_threshold_low ==0.0 : 
            p_threshold_low= 1e-7
            print("correcting its threshold to :",p_threshold_low)

    decoratedFlowRead = gsd.fl.GSDFile(name=filename, mode='rb')   # TODO : use with as instead of close()
    atomTypes = decoratedFlowRead.read_chunk(frame=0, name='atomTypes')
    boxParams   = decoratedFlowRead.read_chunk(frame=0, name='configuration.box')   
    Lx = boxParams[0]
    Ninitial = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'N')))
    expectedDensity = float(module_filenamesHandler.filename_parser(filename[:-4], 'rho'))
    selectedAtoms, Natoms, expectedDensity = module_filenamesHandler.selectAtoms_function(atomTypes, kind, Lx, expectedDensity, 0)
    selectedAtoms = np.arange(Ninitial, dtype=int)[selectedAtoms]
#    selectedAtoms = np.arange(Natoms, dtype=int)
#    if 'tagBonusNum' in filename :
#       Natoms = int(module_filenamesHandler.filename_parser(filename[:-4], 'tagBonusNum')) 
#       selectedAtoms = selectedAtoms[:Natoms]
    atomTypes = atomTypes[selectedAtoms]
    atomIds   = decoratedFlowRead.read_chunk(frame=0, name='atomIds')[selectedAtoms]
    NactMax   = decoratedFlowRead.read_chunk(frame=0, name='NactMax')[0]
    
    
    unique_times_number = decoratedFlowRead.nframes
    unique_times    = np.zeros(unique_times_number, dtype=int)
#    trajectory      = np.zeros((unique_times_number, Natoms, 3))
    phop            = np.zeros((unique_times_number, Natoms))
#    times           = np.zeros((unique_times_number, Natoms), dtype=int)
#    TAs             = np.zeros((unique_times_number, Natoms), dtype=int)
#    TBs             = np.zeros((unique_times_number, Natoms), dtype=int)
#    averageA        = np.zeros((unique_times_number, Natoms, 3))
#    Prev_StdDev     = np.zeros((unique_times_number, Natoms))
#    Next_StdDev     = np.zeros((unique_times_number, Natoms))
#    average_distAB  = np.zeros((unique_times_number, Natoms))
#    softness        = np.zeros((unique_times_number, Natoms))
    for tc_index in range(unique_times_number):
        unique_times[tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='tc')[0]
#        trajectory[tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='positions')[selectedAtoms,:]
        phop     [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='phop')[selectedAtoms]    #  resultsArray[:,tc_index,1]
##        times     [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='times')[selectedAtoms]
#        TAs       [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='TAs')[selectedAtoms]
#        TBs       [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='TBs')[selectedAtoms]
#        averageA  [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='averageA')[selectedAtoms]
#        Prev_StdDev[tc_index]= decoratedFlowRead.read_chunk(frame=tc_index, name='Prev_StdDev')[selectedAtoms]
#        Next_StdDev[tc_index]= decoratedFlowRead.read_chunk(frame=tc_index, name='Next_StdDev')[selectedAtoms]
#        average_distAB[tc_index]= decoratedFlowRead.read_chunk(frame=tc_index, name='average_distAB')[selectedAtoms]
#        softness  [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='softness')[selectedAtoms]
    decoratedFlowRead.close()


    print("... file loaded. Use larger every_recAlways if loading is very slow (you may have stupidly set it to 1).")
    print("NactMax=", NactMax,  "    len(phop)=",len(phop))


#    raise SystemExit

    # selection: excludes the every_recAlways points
    activityState = np.zeros( (Natoms, NactMax, 2), dtype=int)
    '''
    ## 1st colum are tiems
    ## 2nd colulm are label:s -1 (hard) , 0 (unlabelled), 1 (p>p_threshold_low), 2 (p>p_threshold_high) 
    '''
    
    Nact_atom = np.zeros( (Natoms), dtype=int)
    for atom in range(Natoms):
#        mask = (phop[:,atom] > p_threshold_low)
#        interesting_indices = list(np.arange(unique_times_number, dtype=int)[mask])
        nact=0
        for tc_interesting in range(unique_times_number):
            if phop[tc_interesting, atom] > p_threshold_low :
#        for tc_interesting in interesting_indices:
    #            interesting_times = list(unique_times[mask])
                activityState[atom, nact, 0] = tc_interesting
                activityState[atom, nact, 1] = 1
                if phop[tc_interesting, atom] > p_threshold_high :
                    activityState[atom, nact, 1] = 2
                    
    #            activityState[atom, nact, 1] = unique_times[tc_interesting]-unique_times[tc_interesting-1] ## for the first element, it will use [-1], the last element==0 
                nact +=1
        
        ## if no events or only 1(+2 throwBoundaries_soft) event(s) :
        if nact < 2 + 2*throwBoundaries_soft: 
            print( "Atom #",atom," :  only nact=", nact,  " found")
            activityState[atom, nact, 0] = 0
            activityState[atom, nact, 1] = -1
            nact += 1
            activityState[atom, nact, 0] = unique_times_number-1
            activityState[atom, nact, 1] = -1
            nact += 1
        Nact_atom[atom] = nact


    ## automatic computation of the appropriate slowWindow size 
    if slowWindow < 0 :
        timeDiffs_allaccumulated = np.zeros(NactMax*Natoms)    
        timeDiffs_allaccumulated_length = 0
        for atom in range(Natoms):
            Nact = Nact_atom[atom]
            NactiveTimes = Nact +1 
            if NactiveTimes -1 -1 > 0:      # - 2*throwBoundaries >0 :
                timeDiffs = np.diff(activityState[atom, :Nact, 0]) # [throwBoundaries:-throwBoundaries]
                timeDiffs_allaccumulated[timeDiffs_allaccumulated_length:timeDiffs_allaccumulated_length+NactiveTimes-1-1] = timeDiffs #  +throwBoundaries: -1-1-throwBoundaries
                timeDiffs_allaccumulated_length += NactiveTimes-1 -1 #- 2*throwBoundaries
        timeDiffs_allaccumulated = timeDiffs_allaccumulated[0:timeDiffs_allaccumulated_length]

        timeDiffs_allaccumulated = np.sort(timeDiffs_allaccumulated)
        ## we take the last Nsamples elements of this sorted array 
        ## to be the large windows. 
        ## The smallest value is also the thresold window size value
        slowWindow = int(timeDiffs_allaccumulated[-Nsamples-1])
        print("\nauto-select slowWindow= ", slowWindow, "  ,  so as to have at least   Nsamples=", Nsamples,"   of independent 'slow' windows (and thus that number of soft instants). \nREMARK: there is no gurantee about the number of soft particles you will find... you need to try. \n")

    ## automatic (estimate) of the correct p_threshold_high (a higher bound)
    if p_threshold_high < 0 :
        p_threshold_high        = np.round(np.sort(phop[phop>0])[-Nsamples    -1],3)
        p_threshold_high_low    = np.round(np.sort(phop[phop>0])[-Nsamples*1.5-1],3)
        print("\nauto-select p_threshold_high< ", str(p_threshold_high), "  maybe down to ", str(p_threshold_high_low)," ,  so as to have at least   Nsamples=", Nsamples,"   instants with p>p_threshold_high. But this is without considering the fastWindow>0, so the p_threshold_high you need may be smaller. (it is exact only in the case fastWindow=0).  \nREMARK: there is no gurantee about the number of hard particles you will find... you need to try. \n")

    rootname_AVA = filename[:-4]+"_kind="+kind+"_pl="+str(p_threshold_low)+"_ph="+str(p_threshold_high)+"_sW="+str(slowWindow)+"_fW="+str(fastWindow)+"_AVA="
    ## TODO : just wriote type=trajectory  and not  -decorated. 
    
    ###################################
    if plotting == True  and slowWindow < 0 :
        plt.figure(40,[20,6])
        base = 1.1
        low_val = 15
        Log_bins_integers = np.array( range(low_val) + [low_val*base**i for i in range(0, int(np.log(np.max(timeDiffs_allaccumulated)/low_val)/np.log(base)), 1) ])-0.000001
        heights, trash = np.histogram(timeDiffs_allaccumulated, bins=Log_bins_integers)
        heights = heights / np.diff(Log_bins_integers)    ## adjust density by histo bin width
        plt.loglog(Log_bins_integers[:-1], heights, lw=3, label='$p_th='+str(p_threshold_low)+'$')
        plt.xlabel(r'$\Delta t = t_{hop,2}-t_{hop,1}$')
        plt.ylabel(r'$N( \Delta t )$')
        plt.legend(loc='best')
        outName=rootname_AVA + "distro-interhopTimes"
        module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
        if debug_or_choose_parameters == False:
            plt.close(40)
    ###################################

    
    length_array_record = int(1e7)
    hard_windows = np.zeros( (length_array_record, 3), dtype=int)
    hard_num = 0
#    med_windows  = np.zeros( (length_array_record, 3), dtype=int)
#    med_num = 0
    soft_windows = np.zeros( (length_array_record, 3), dtype=int)
    soft_num = 0

    ## hard windows:  when p<p_threshold_low for a time slowWindow or more  ##
    for atom in range(Natoms):
        for nact in range(1,Nact_atom[atom]):
            ## activityState is defined only when phop > p_low
            if activityState[atom, nact, 0] - activityState[atom, nact-1, 0] >= slowWindow : 
                hard_windows[hard_num, 0] = atom
                hard_windows[hard_num, 1] = activityState[atom, nact-1, 0]
                hard_windows[hard_num, 2] = activityState[atom, nact  , 0] 
                hard_num += 1
    hard_windows = hard_windows[:hard_num]      
               
               
    ## TODO: jit this ##
    ## soft windows: when p>p_threshold_high for some time, or below it only very shorty (less than fastWindow of time) ##
    for atom in range(Natoms):
        nact = 0  
        while nact < Nact_atom[atom]-1 : 
            t_prev = activityState[atom, nact  , 0]
            reallyActive = 0
            nact_initial = nact
            if activityState[atom, nact, 1] == 2:
                    reallyActive = 1
            while activityState[atom, nact, 0]+1 == activityState[atom, nact+1, 0] and nact+1 < Nact_atom[atom]-1  :
                if activityState[atom, nact, 1] == 2:
                    reallyActive = 1
                nact+=1
            if reallyActive == 1:
                soft_windows[soft_num, 0] = atom
                soft_windows[soft_num, 1] = activityState[atom, nact_initial, 0] # begining of soft window 
                soft_windows[soft_num, 2] = activityState[atom, nact, 0]
                soft_num += 1
            nact += 1

#    for atom in range(Natoms):
#        nact = 0   
#        nact = 0
#        while nact < Nact_atom[atom]:
#            if activityState[atom, nact, 1] == 2 :
#                soft_windows[soft_num, 0] = atom
#                soft_windows[soft_num, 1] = activityState[atom, nact, 0]# begining
#                active = 1
#                while active :
#                    while activityState[atom, nact, 1] == 2:
#                        nact += 1
#                    nact -= 1

#                    if activityState[atom, nact+1, 0] - activityState[atom, nact, 0] < fastWindow : 
#                        if activityState[atom, nact+1, 1] == 2 :
#                            nact += 1 # and we stay at active=1
#                        else :
#                            active = 0
#                    else:
#                        active = 0
#                        
#                soft_windows[soft_num, 2] = activityState[atom, nact, 0]    
#                soft_num += 1
#                ## end of the soft window ##
#            nact +=1
    soft_windows = soft_windows[:soft_num]

    


################################################################################



#    #### DEBUG ####
#    atom=45
#    NactiveTimes = activeTimes[atom,0]
#    timeDiffs = np.diff(activeTimes[atom,1:NactiveTimes])
#    print(AtomTimesHard[AtomTimesHard[:,0]==atom])
#    print(AtomTimesSoft[AtomTimesSoft[:,0]==atom])
#    print("timeDiffs",timeDiffs)
#    print("activeTimes[atom,:NactiveTimes+2]",activeTimes[atom,:NactiveTimes+2])
#    plt.figure(501)
#    plt.plot(activeTimes[atom,1:NactiveTimes-1], timeDiffs) 
#    ## spike of this are the long quiescent times 
#    #### DEBUG ####


    print("we found hard_num=", hard_num, "  hard particles.")
    print("we found soft_num=", soft_num, "  soft particles.")
    #check there are enough big events (at least as many as Nsoft used )

    if soft_num < Nsamples or hard_num < Nsamples:
        if hard_num < Nsamples :
            print("Warning: Not enough HARD space-time points found: only ", hard_num) 
            print("Play with the option slowWindow==-1 to find a decent compromise between Nsamples and slowWindow")
        if soft_num < Nsamples :
            print("Warning: Not enough SOFT space-time points found: only ", soft_num) 
            print("Play with the option p_threshold_high to find a decent compromise between Nsamples and p_threshold_high")
        print("\nYou have NOT found your favorite slowWindow+p_threshold_low/p_threshold_high+fastWindow/Nsamples compromise!  This should be because of not enough soft particles, as the slowWindow was automatically computed so as to have enough hard particles.\nRe-run this program using a smaller  p_threshold_high<",p_threshold_high," \n\n")
#        raise SystemExit
#        min(hard_num, soft_num)
    else: 
        ## TODO: simplify filenames, by using pl and slowWindow for the hard particles 
        ##                                      and only ph and fastWindow  ofr the soft.
        header= "Displays the instants (atom+time) of the beggining (begg) and end of each window that is identified as soft/hard (see the file name)\natomIds  begg  end"
        
        softwindowsName = rootname_AVA + "trainers-soft.dat"
        np.savetxt(softwindowsName, soft_windows, fmt="%d %d %d", header=header)
        print(softwindowsName, "    was written to a file")

        hardwindowsName = rootname_AVA + "trainers-hard.dat"
        np.savetxt(hardwindowsName, hard_windows, fmt="%d %d %d", header=header)
        print(hardwindowsName, "    was written to a file")
        
################################################################################################################################################################


    if 'SLID' not in filename:
        print("\n CAREFUL !! jai pas encore adate cet algo a la version IBS de phop .  Donc les fenetres rapides sont pas forcement ok? Ou bien les seuils phop sont a ajuster")










################################################################################################################################################################
## TODO: delete all this  crap:
#    ## record the universe instants (space==atom number ; time step) that are interesting ##
#    AtomTimesHard = np.zeros((length_array_record, 3), dtype=int)
#    hard_num = 0
#    AtomTimesSoft = np.zeros((length_array_record, 3), dtype=int)
#    soft_num = 0


#    for atom in range(Natoms):   
#        ## each atom has a different number of events ##
#        NactiveTimes = int(resultsArray[atom, 0 ,11])
#        activeTimes = np.array(resultsArray[atom,0:NactiveTimes,11],dtype=int)
##        resultsArray[atom,1:NactiveTimes,11] = np.sort(resultsArray[atom,1:NactiveTimes,11])        


#        activityState = np.zeros( (NactiveTimes,2), dtype=int )
#        activityState[:,0] = activeTimes
#        
#        activityState[:,0] = 0 ## the  first time  is artificially set to be 0
#        activityState[:,1] = 1 ## that first event is artificially set to be active (1)
#        t0 = 0
#        j = 0
#        for i in range(1,NactiveTimes) : 
#            t1 = activeTimes[i]
#            if t1-t0 > delay_soft  : 
#                activityState[j,0] = t0
#                activityState[j,1] = 0
#                activityState[j+1,0] = t1
#                activityState[j+1,1] = 1
#                j += 2
#            t0 = t1
#        activityState = activityState[:j,:]
#        f.write_chunk(name='activityState',  data=activityState ) # in untis of "tc"


#####        douteux: 
#    Delta_r        = np.zeros((Natoms, NactMax))
#    phop_recomputed= np.zeros((Natoms, NactMax))
#    Var_previous   = np.zeros((Natoms, NactMax))
######        Delta_r        [atom,0:NactiveTimes] = eventsDataFlowRead.read_chunk(frame=atom, name='Delta_r')
######        phop_recomputed[atom,0:NactiveTimes] = eventsDataFlowRead.read_chunk(frame=atom, name='phop_recomputed')
######        Var_previous   [atom,0:NactiveTimes] = eventsDataFlowRead.read_chunk(frame=atom, name='Var_previous')


#    ## TODO: do I need gsd for this ? 
#    ## maybe recording scipy sparse matries would be enguh ?
#    gsd.fl.create(name=eventsDataName,application="events Data, ie atom per atom file",\
#            schema="each frame is an atom. Data chunks are undesrtood over time",\
#            schema_version=[1,0])
#    f = gsd.fl.GSDFile(name=eventsDataName, mode='wb')
#    # writing atom types only in the first frame : (and not closing this first frame) #
#    f.write_chunk(name='atomTypes' , data=np.array(atomTypes , dtype=int))
#    f.write_chunk(name='configuration.box' , data=np.array(boxParams , dtype=np.float32))
#    f.write_chunk(name='atomIds' , data=np.array(atomIds , dtype=int))
#    f.write_chunk(name='curStep'  , data=np.array([curStep], dtype=int) )
#    for atom in range(Natoms):  
#        f.write_chunk(name='activityState',    data=activityState   ) # in untis of "tc"
#        f.end_frame()
#    f.close()
#    print("finished writing  ", eventsDataName,  " to file.")

