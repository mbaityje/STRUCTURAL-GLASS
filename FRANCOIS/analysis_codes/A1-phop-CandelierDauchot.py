from __future__ import print_function
import gsd.pygsd
import gsd.hoomd
import gsd.fl

import sys
import numpy as np

import module_filenamesHandler 
import module_phop


######
## plots: 
import module_importPlotParams
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True          
plt.rcParams.update(module_importPlotParams.fparams())
plt.ion()   ## optional
######



# TODO: compare :
# hops found between this, and when you re-compute p_hop using the tc of the neighboring p_hops as TA,TB boundaries.
## stuffs may become more homoegenous, like the very large phop would be out, and then 
## utiliser le tc-1 et tc+1 de chaque tc donne en tant que valeurs TA Tb, pour voir si ca change le phop trouve au milieu ?

# TODO: voir est ce que IBS mets plus de sauts pour les particules B ? Les 200 dernieres?
########################

# TODO: run this code on realistic input with higher density.
# analyze the phop distro output, select a decent phop, 
# do a full study on p_threshold_low dependence

# TODO:  plots of softness. .. everything bascially (assuming softness is recorded there)


print('usage:')
print('run  phop*.py    [filename]    [frameSkipEvery=1]   [p_threshold=0.002] ')
filename = sys.argv[1]
rootname = module_filenamesHandler.get_rootname(filename)
suffix   = module_filenamesHandler.get_suffix(filename)

frameSkipEvery=1
if len(sys.argv)>2:
    frameSkipEvery = int(sys.argv[2])

p_threshold_low  = 0.002
if len(sys.argv)>3:
    p_threshold_low = float(sys.argv[3])

p_threshold_high = 1.0
every_recAlways = 100    ## TODO: repair the plots viz. this 
plotting=1

delay_soft=3


slidingWindow=False
slidingWindow=True
IBS = True - slidingWindow


#while(it_I != States_I.end())
#{
#        var_I += (*it_I - MeanPositionJ).segment<DIM>(DIM*i).squaredNorm()/States_I.size();
#        var_J += (*it_J - MeanPositionI).segment<DIM>(DIM*i).squaredNorm()/States_J.size();                                
#        //if(i==1177)
#        //cout << boxsize*(*it_I).segment<DIM>(DIM*i).transpose() << "\n";
#        it_I++;
#        it_J++;
#}
#double p_hop = boxsize*boxsize*sqrt(var_I*var_J);
## boxsize == Lx -> convert to real units
## 

if "deco" not in filename and ('traj' in filename or 'FIRE' in filename): 

################################################################################
####### read the file (using the gsd pakage) ##################################
    with open(filename, 'rb') as flow: 
        HoomdFlow= gsd.pygsd.GSDFile(flow)
        hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
        dim         = hoomdTraj.read_frame(0).configuration.dimensions    
        Nframes     = hoomdTraj.file.nframes
        boxParams   = hoomdTraj.read_frame(0).configuration.box
        Natoms      = hoomdTraj[0].particles.N
        atomIds     = np.arange(Natoms, dtype=int).transpose() ## XXX to be added ?
        atomTypes   = np.array(hoomdTraj[0].particles.typeid, dtype=int).transpose()
        Lx=boxParams[0]
        if boxParams[0] != boxParams[1] or  boxParams[1]!= boxParams[2] or  boxParams[0]!= boxParams[2]:
            print( 'box dimensions are : ', boxParams[0])
            print( 'and are not supported (only isotropic systems supported)')
            raise SystemExit

        ## TODO: use this time unit in the outputs ??
        dt=0.00025
        recPeriod = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'rP')))
        curStep = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'cst')))
        time_per_index = recPeriod*dt * frameSkipEvery # time_per_frame = recPeriod * dt
        
        Nframes = len(hoomdTraj)/frameSkipEvery
        Nframes_max = int(6e8/Natoms/3)
        trajDuration = min(Nframes_max, Nframes)
        ## trajDuration should be related to the actual time observed in the simu, 
        ## or at least the number of steps seen, indpeenent of number of atoms in the simulation...
        print('there are Nframes=', Nframes, 'in the file, but we only use trajDuration=', trajDuration, ' of them.')

    #    print('TODO: if you need longer time sequences, you can split the atoms instead of splitting time, and do separate chuncks of the system using a mask, like selectedAtoms')
        ## for tests and/or RAM limitations:
        selectedAtoms = np.arange(Natoms, dtype=int) #[:50] 
        additional_tag='-'
        if slidingWindow:
            additional_tag += "SLID"
        elif IBS:
            additional_tag += "IBS"
#        if len(selectedAtoms) < Natoms:
#            additional_tag += "-tgBnus="+str(len(selectedAtoms))
        
        Natoms = len(selectedAtoms)
        trajectory = np.zeros((trajDuration,Natoms,3))
        for time in range(0, trajDuration, 1):
            trajectory[time] = (hoomdTraj[time*frameSkipEvery].particles.position)[selectedAtoms]
        atomIds = atomIds[selectedAtoms]
        HoomdFlow.close() # is closed automatcally here ##
##### DEBUG: this allows to check windows B and A are equivalently treated, exactly : 
#        trajectory = trajectory[::-1,:]
#        print 'inverted order'
##### end of debug ###
        
    print('Shape of your trajectory array:', np.shape(trajectory))
    ## "trajectory" is now the full array of all positions of all atoms (over time)
    ## trajectory[:,1] is the trajectory of particle 1 over time [its  (x,y,z) coords t be precise]
################################################################################





    print('searching for hops larger than p_thres')

    if slidingWindow == True:
        resultsArray, Nact_atom, hop_number = module_phop.FixedWindow_phopFromTrajectory_func(trajectory, Lx, p_threshold_low, every_recAlways)
        
    elif IBS == True:
        resultsArray, Nact_atom, OldIntervals_lengthUsed_max, hop_number = module_phop.phopFromTrajectory_func(trajectory, Lx, p_threshold_low, every_recAlways)
        print("the maximal width of the tree was ", OldIntervals_lengthUsed_max, "  (increase OldIntervals_lengthMax if you have memory errors.. before this point)")

    print('finished searching for events. #events found=',hop_number)



################################################################################
    ## get rid of multiplicity, and of the frames recorded every_recAlways ##
    unique_times = resultsArray[(resultsArray[:,:,1] > p_threshold_low), 0]
    # destroys all duplicates, using a hash method: #
    unique_times = np.sort(np.array(list(set(list(unique_times))), dtype=int))
    unique_times_number = len(unique_times)
    print('we found ', unique_times_number, ' unique time of phops    with p >', p_threshold_low)


    ## TODO: actually compute softness here, on-the-fly
    softness = np.zeros(Natoms,  dtype=np.float32)

    ## creer le fichier source avec le new name: 
    trajType = str(module_filenamesHandler.filename_parser(filename[:-4], 'type'))
    decoratedName =rootname+"_type=" +trajType+'-deco'+additional_tag+"_"+suffix[:-4]+"_pLow="+str(p_threshold_low)+"_rec="+str(every_recAlways)+""
    if frameSkipEvery != 1:
        decoratedName += "_skp="+str(frameSkipEvery)
    decoratedName += ".gsd"
#########    eventsDataName=rootname+"_type="+trajType+'-eventsData_'+suffix[:-4]+additional_tag+"_pLow="+str(p_threshold_low)+"_skip="+str(frameSkipEvery)+"_rec="+str(every_recAlways)+""+".gsd"
    

    
    gsd.fl.create(name=decoratedName,application="trajectory with decorations",\
            schema="[atomTypes] [atomIds] tc positions phop TAs TBs averageA Prev_StdDev Next_StdDev average_distAB softness",\
            schema_version=[1,0])
    f = gsd.fl.GSDFile(name=decoratedName, mode='wb')
    # writing atom types only in the first frame : (and not closing this first frame) #
    f.write_chunk(name='atomTypes' , data=np.array(atomTypes , dtype=int))
    f.write_chunk(name='configuration.box' , data=np.array(boxParams , dtype=np.float32))
    f.write_chunk(name='atomIds' , data=np.array(atomIds , dtype=int))
    f.write_chunk(name='curStep'  , data=np.array([curStep], dtype=int) )
    
    NactMax = np.max(Nact_atom) # resultsArray[:, 0 ,11])
    f.write_chunk(name='NactMax'  , data=np.array([NactMax], dtype=int) )
    
    for tc_index in range(unique_times_number):
        tc = unique_times[tc_index]
        f.write_chunk(name='tc'  , data=np.array([tc], dtype=np.float32) ) ## DEBUG ONLY
        # we rebuild the positions in a separate file, becasue not all frames are present
        f.write_chunk(name='positions' , data=np.array(trajectory[tc], dtype=np.float32))
        f.write_chunk(name='softness' , data=softness )

        ### to be recorded... in particular after second sweep
        f.write_chunk(name='phop'     , data=np.array(resultsArray[:,tc,1]  , dtype=np.float32))
        
        
#        ## TODO : compute the Delta_r relative to each hop, put it here, using the mid-points *instanteeous* positions,
#        f.write_chunk(name='Delta_r',         data=np.array( [1]  , dtype=np.float32))  ## resultsArray[atom,activeTimes, 1]
## RQ XXX : set Deltar to 0 for the atoms that did not move ! Since it's deltaR in between windows (across one hop), not within a dt or 1 frame. 
#        f.write_chunk(name='phop_recomputed', data=np.array( [1]  , dtype=np.float32)) ## TODO : compute it
    ## TODO: use the knowledge of Delta_r or of phop_recomputed in the next script, construct_TrainingSet.py ## 

 
#TODO :  remove element 2 ??
#        f.write_chunk(name='event'    , data=np.array(resultsArray[:,tc,2]  , dtype=int))## bool is not supported by the gsd record types
        
    ## the following elements are mostly needed for phop-debugging. In the end we do not need them...
    ## we want to know values AFTER the second sweep, not before (they are unphysical)
        ## note that TAs and TBs  correspond only to the TA TB of the phop when it was recorded, but further events may have occured closer than this in time.
        f.write_chunk(name='TAs', data=np.array(resultsArray[:,tc,3]  , dtype=int))
        f.write_chunk(name='TBs', data=np.array(resultsArray[:,tc,4]  , dtype=int))
        averageA = np.array([resultsArray[:,tc,5], resultsArray[:,tc,6], resultsArray[:,tc,7] ]  , dtype=np.float32).transpose()
        f.write_chunk(name='averageA'   , data=np.array(averageA              , dtype=np.float32))
        f.write_chunk(name='Prev_StdDev', data=np.array(resultsArray[:,tc,8]  , dtype=np.float32))
        f.write_chunk(name='Next_StdDev', data=np.array(resultsArray[:,tc,9]  , dtype=np.float32))
        f.write_chunk(name='average_distAB', data=np.array(resultsArray[:,tc,10]  , dtype=np.float32))
        f.end_frame();
    print("number of frames recorded to a file: ", f.nframes)
    f.close()
    #f.write_chunk(name='xi', data=np.array(resultsArray[:,tc,2], dtype=np.float32))
    #xi = decoratedFlowRead.read_chunk(0,'xi')
    print("finished writing  ", decoratedName,  " to file.\n\n")

    
    ## now, follow up with the analysis ! ##
    filename = decoratedName
################################################################################




################################################################################
#### reading the "decorated" file, and plot stuff #####
if "deco" in filename :   # read the decorated file produced before # 

    p_threshold_low = float(module_filenamesHandler.filename_parser(filename[:-4], 'pLow'))
    print('reading that the file has a p_threshold_low=', p_threshold_low, ' ...')
    if p_threshold_low ==0 : 
        p_threshold_low= 1e-7
        print("correcting its threshold to :",p_threshold_low)

    Natoms = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'N')) )
    selectedAtoms = np.arange(Natoms, dtype=int)
    if 'tagBonusNum' in filename :
       Natoms = int(module_filenamesHandler.filename_parser(filename[:-4], 'tagBonusNum')) 
       selectedAtoms = selectedAtoms[:Natoms]

    decoratedFlowRead = gsd.fl.GSDFile(name=filename, mode='rb')   # TODO : use with as instead of close()
    atomTypes = decoratedFlowRead.read_chunk(frame=0, name='atomTypes')[selectedAtoms]
    atomIds   = decoratedFlowRead.read_chunk(frame=0, name='atomIds')[selectedAtoms]
    boxParams   = decoratedFlowRead.read_chunk(frame=0, name='configuration.box')   
    Lx = boxParams[0]
    Natoms = len(atomTypes)
    unique_times_number = decoratedFlowRead.nframes
    unique_times    = np.zeros(unique_times_number, dtype=int)
    trajectory      = np.zeros((unique_times_number, Natoms, 3))
    phops           = np.zeros((unique_times_number, Natoms))
    times           = np.zeros((unique_times_number, Natoms), dtype=int)
    TAs             = np.zeros((unique_times_number, Natoms), dtype=int)
    TBs             = np.zeros((unique_times_number, Natoms), dtype=int)
    averageA        = np.zeros((unique_times_number, Natoms, 3))
    Prev_StdDev     = np.zeros((unique_times_number, Natoms))
    Next_StdDev     = np.zeros((unique_times_number, Natoms))
    average_distAB  = np.zeros((unique_times_number, Natoms))
##    softness        = np.zeros((unique_times_number, Natoms))
    for tc_index in range(unique_times_number):
        unique_times[tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='tc')[0]
        trajectory[tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='positions')[selectedAtoms,:]
        phops     [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='phop')[selectedAtoms]    #  resultsArray[:,tc_index,1]
##        times     [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='times')[selectedAtoms]
        TAs       [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='TAs')[selectedAtoms]
        TBs       [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='TBs')[selectedAtoms]
        averageA  [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='averageA')[selectedAtoms]
        Prev_StdDev[tc_index]= decoratedFlowRead.read_chunk(frame=tc_index, name='Prev_StdDev')[selectedAtoms]
        Next_StdDev[tc_index]= decoratedFlowRead.read_chunk(frame=tc_index, name='Next_StdDev')[selectedAtoms]
        average_distAB[tc_index]= decoratedFlowRead.read_chunk(frame=tc_index, name='average_distAB')[selectedAtoms]
##        softness  [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='softness')[selectedAtoms]
    decoratedFlowRead.close()
    rootname_PHOP=filename[:-4]+"_PHOP="


    print("... file loaded. Use larger every_recAlways if loading is very slow (you may have stupidly set it to 1).")

    # selection: excludes the every_recAlways points
    basic_mask = (phops>p_threshold_low)  ## i.e. all of them

    every_plot=100
    
    Ttots = TBs-TAs

    Cross_StdDev = Prev_StdDev * Next_StdDev
    Cross_StdDev[Cross_StdDev==0] +=0.1  # TODO : il  y a des divergences lorsque trop de sauts on tpas ete vus ou bien ont une vraiance de 0 (proabbement a casue du fait qu el'intervalle est vide !!)

    tstoplot = unique_times#[basic_mask]    # all the times (non unique with events) 
    decent_phop = phops[basic_mask]
#    


#    ########## phop VS time  ############

#    ## phop(t) for a few atoms ##
#    plt.figure(11,[10,6])
#    every_plot=1
#    kind = 1 ## plot the p_hop_value
#    ts = np.arange( 0, unique_times_number,every_plot ,dtype=int)
#    for atom in range(45,50):
#        plt.semilogy(ts, phops[::every_plot, atom], ls='-', marker='')
#        plt.semilogy(ts, Cross_StdDev[::every_plot, atom], ls='-.', marker='')
#    plt.xlabel(r'$t$')
#    plt.ylabel(r'$p_{hop}$')
#    outName=rootname_PHOP + "several-atoms-phop"
#    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())

    ## plot the p_hop_value's   for a single atom ##
    plt.figure(10,[10,6])
    every_plot=1
    atom = 44
#    ts = np.arange( 0, unique_times_number ,dtype=int)
#    mask = phops[:,atom]>0
    mask_nonzero = (phops[:,atom]>p_threshold_low)
    hoptoplot = np.log(10.0*phops[mask_nonzero,atom])
    plt.plot((unique_times[mask_nonzero])[::every_plot], hoptoplot[::every_plot], ls='', marker='o', lw=2.0)
    unrolled_traj  = module_phop.unrolled_traj_func(trajectory, Lx, atom)
    unrolled_traj -= np.mean(unrolled_traj,0)
    plt.plot(unique_times[::every_plot],unrolled_traj[::every_plot]*3.0)
    plt.title(r'atom $'+str(atom)+'$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$3x,3y,3z, \log(p_{hop})$')
    outName=rootname_PHOP + "single-atom-traj-"+str(atom)
    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
    plt.ylim([-4,2])
    outName=rootname_PHOP + "several-atoms-phop-cut"
    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())

#plt.plot(activityState[:,0], activityState[:,1])       ## XXX 

#    ### phops of all atoms over time: ####
#    ### (scatterplot of all activity)
#    ## selecting only the big phop events
#    ## we flatten the array at this point (Natoms+trajDuration)
#    plt.figure(13,[10,6])
#    print('we found ', len(decent_phop), ' decent_phop with   p >', p_threshold_low)
### scatterplot: decent_phop(t) for the decent_phop, all atoms.
#    plt.semilogy(tstoplot, decent_phop, ls='', marker='x')
#    plt.xlabel(r'$t_{event}$')
#    plt.ylabel(r'$p_{hop}$')
##    outName=rootname_PHOP + "scatterplot_all" #__phop_geq"+str(p_threshold_low)
#    outName=rootname_PHOP + "scatterplot_phop_VS_time" #_geq"+str(p_threshold_low)
#    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())

    ### (density of activity over time)
    ## plotting the times of events (and how many events there are at each of these times)
    plt.figure(30,[10,6])
    activity, bins = np.histogram(tstoplot, bins = int(unique_times_number/10.0)+2)
    plt.xlabel("$t_c$ (in the frames+frameSkipEvery time units) ")
    plt.ylabel("#hops$(t_c)$")
    plt.plot(bins[:-1], activity, ls='-', marker='')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$N_{event}$ (activity)')
    outName=rootname_PHOP + "activity_VS_time" #_geq"+str(p_threshold_low)
    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
    ######################################################



    ###### distribution of p_hop_value ######
    base = 1.1
    Log_bins = np.array([p_threshold_low*base**i for i in range(-5, int(np.log(np.max(decent_phop)/p_threshold_low)/np.log(base)), 1) ])
    heights, trash = np.histogram(decent_phop, bins=Log_bins)
    heights = heights / np.diff(Log_bins)    ## adjust density by histo bin width
    plt.figure(40,[20,6])
    plt.loglog(Log_bins[:-1], heights, lw=3)
    plt.xlabel(r'$p_{hop}$')
    plt.ylabel(r'$N(p_{hop} | p_{hop}>'+str(p_threshold_low)+')$')
    outName=rootname_PHOP + "distro-phop-loglog_phop_geq"+str(p_threshold_low)
    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
#    a.append((Log_bins[:-1], heights))

#    ## same in semilogy 
#    plt.figure(41,[10,6])
#    plt.semilogy(Log_bins[:-1], heights, lw=3)
#    plt.xlabel(r'$p_{hop}$')
#    plt.ylabel(r'$N(p_{hop} | p_{hop}>'+str(p_threshold_low)+')$')
#    outName=rootname_PHOP + "distro-phop-semilogy_phop_geq"+str(p_threshold_low)
#    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())




#    ################ phop VS Ttot ################

##    ## scatterplot
##    plt.figure(15,[10,6])
##    plt.loglog(phops[basic_mask], Ttots[basic_mask], ls='', marker='x')
##    plt.xlabel(r'$p_{hop}$')
##    plt.ylabel(r'$T_f$')
##    outName=rootname_PHOP + "phop_Ttot"
##    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())

##    ## scatterplot
##    plt.figure(16,[10,6])
##    plt.loglog(Ttots[basic_mask], phops[basic_mask], ls='', marker='x')
##    plt.xlabel(r'$T_f$')
##    plt.ylabel(r'$p_{hop}$')
##    outName=rootname_PHOP + "Ttot_phop"
##    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
#    
#    
    plt.figure(1,[20,6])
    plt.figure(2,[20,6])
    plot_histograms = 1
    if plot_histograms:
        ## P(p_{hop} | T_{tot})
        averagePhop=[]
        base_Tflist=4.0
        Tflist=[3,10,30,100,300,1000,3000]
        Tflist=[2]+[base_Tflist**i for i in range(1,7)]
        base=1.1
        Log_bins = np.array([p_threshold_low*base**i for i in range(-40, int(np.log(np.max(phops)/p_threshold_low)/np.log(base)), 1) ])
        for Tf in Tflist:
            mask = ((Ttots > Tf) * (Ttots < Tf*base_Tflist) * (phops>p_threshold_low))
            phops_tobin = phops[mask]
            Nitem = len(phops_tobin)
            averagePhop.append(np.mean(phops_tobin))
            if Nitem>1: 
                P_phop, trash = np.histogram(phops_tobin, bins=Log_bins)
                P_phop = P_phop / np.diff(Log_bins)   ## adjust density by histo bin width
                plt.figure(1,[20,6])
                plt.loglog(Log_bins[:-1], P_phop        , lw=3, label=r'$T_{tot}\in ['+str(Tf)+', '+str(Tf*base_Tflist)+']$')
                plt.figure(2,[20,6])
                plt.loglog(Log_bins[:-1], P_phop / Nitem, lw=3, label=r'$T_{tot}\in ['+str(Tf)+', '+str(Tf*base_Tflist)+']$')
            
        plt.figure(1,[20,6])
        plt.legend(loc='best')
        plt.xlabel(r'$p_{hop}$')
        plt.ylabel(r'$N(p_{hop} | T_{tot}\in [ ... ] )$')
        outName=rootname_PHOP + "distroN-phop_of_Ttot"
        module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())

        plt.figure(2,[20,6])
        plt.legend(loc='best')
        plt.xlabel(r'$p_{hop}$')
        plt.ylabel(r'$P(p_{hop} | T_{tot}\in [ ... ] )$')
        outName=rootname_PHOP + "distroP-phop_of_Ttot"
        module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())


        ## < (p_{hop} | T_{tot}) >
        plt.figure(18,[10,6])
        plt.loglog(Tflist, averagePhop, marker='o')
        plt.legend(loc= 'best')
        plt.xlabel(r'$T_{tot}$')
        plt.ylabel(r'$<p_{hop}>(T_{tot})$')
        outName=rootname_PHOP + "average-phop-of-Ttot"
        module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
##    ## distribution of p_hop_value (for any Ttot value > 5 == Ttot_threshold) ##
##    plt.figure(40,[10,6])
##    Ttot_threshold=5
##    mask = (Ttots > Ttot_threshold)
##    phops_tobin = phops[mask]
##    
##    base=1.1
##    Log_bins = np.array([p_threshold_low*base**i for i in range(-40, int(np.log(np.max(phops_tobin)/p_threshold_low)/np.log(base)), 1) ])
##    
##    P_phop, trash = np.histogram(phops_tobin, bins=Log_bins)
##    P_phop = P_phop / np.diff(Log_bins)    ## adjust density by histo bin width
##    plt.loglog(Log_bins[:-1], P_phop, lw=3)
##    plt.xlabel(r'$p_{hop}$')
##    plt.ylabel(r'$N(p_{hop} | T_{tot}>'+str(Ttot_threshold)+')$')
##    outName=rootname_PHOP + "distro-phop--forall_Ttot_geq"+str(Ttot_threshold)
##    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())



#    ########## phop VS variances ############

    figNums=70
    XXarray = Cross_StdDev
    Xstr='\sigma'
    XstrNoLatex='sigma'
    base = 1.3
    val_typ = 0.0000001
    Log_bins = np.array([1e-9]+[val_typ*base**i for i in range(-5, int(np.log(np.max(phops)/val_typ)/np.log(base)), 1) ])

#    plt.figure(figNums+0,[10,6])
##    plt.loglog(Cross_StdDev[::10], phops[:-1:10], ls='', marker='x')
#    plt.loglog(Cross_StdDev[basic_mask], phops[basic_mask], ls='', marker='x')
#    plt.xlabel(r'Std. Dev. Cross.')
#    plt.ylabel(r'$p_{hop}$')
#    outName=rootname_PHOP + "scatterplot_phop_VS_StdDev"
#    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
#    plt.close()

    
    plot_histograms = 1
    if plot_histograms:
        base_XXflist=10
        XXflist=[1e-8*base_XXflist**i for i in range(0,9)]
        ## P(p_{hop} | XX_{tot})
        averagePhop=[]
        plt.figure(figNums+1,[20,6])
        base=1.1
        for XXf_index in range(len(XXflist)-1):
            XXf = XXflist[XXf_index]
            mask = ((XXarray > XXf) * (XXarray < XXflist[XXf_index+1]))
            phops_tobin = phops[mask]
            Nitem = len(phops_tobin)
            averagePhop.append(np.mean(phops_tobin))
            if Nitem >1 :
                P_phop, trash = np.histogram(phops_tobin, bins=Log_bins)
                P_phop = P_phop / np.diff(Log_bins) / Nitem   ## adjust density by histo bin width
                plt.loglog(Log_bins[:-1], P_phop, lw=3, label=r'$'+Xstr+'_{tot}\in ['+str(XXf)+', '+str(XXflist[XXf_index+1])+']$')
        plt.legend(loc='best')
        plt.xlabel(r'$p_{hop}$')
        plt.ylabel(r'$N(p_{hop} | '+Xstr+'_{tot}\in [ ... ] )$')
        outName=rootname_PHOP + "distro-phop_of_"+XstrNoLatex
        module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())


## XXX sometimes empty stuff makes this wrong !!
        ## < (p_{hop} | XX_{tot}) >
        plt.figure(figNums+2,[10,6])
        plt.loglog(np.array(XXflist[:-1])*base_XXflist**0.5, averagePhop, marker='o')
        plt.legend(loc= 'best')
        plt.xlabel(r'$'+Xstr+'_{tot}$')
        plt.ylabel(r'$<p_{hop}>('+Xstr+'_{tot})$')
        outName=rootname_PHOP + "average-phop-of-"+XstrNoLatex
        module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())


    heights, trash = np.histogram(Cross_StdDev, bins=Log_bins)
    heights = heights / np.diff(Log_bins) / len(Cross_StdDev)    ## adjust density by histo bin width
    plt.figure(figNums+4,[20,6])
    plt.loglog(Log_bins[:-1], heights, lw=3)
    plt.xlabel(r'$'+Xstr+'$')
    plt.ylabel(r'$P('+Xstr+')$')
    outName=rootname_PHOP + "distro-"+XstrNoLatex
    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())



########### unavailable data (because numba does not handle lists... and I did not coeded for this unessential, algorithmic-dependent feature)
######    ## now plotting histogram of each depth of the phop search: we expected larger phops for initial searches and reducing phop value as you go down the ladder, but...
######    #plt.figure(1111,[10,6])
######    depth_max = 10
######    for p_hop_hierarchy_depth in range(depth_max):
######        heights, trash = np.histogram(p_hop_hierarchy[p_hop_hierarchy_depth],Log_bins)
######        plt.loglog(Log_bins[:-1],heights, lw=(depth_max-p_hop_hierarchy_depth)/2.0)

################################################################################
















#################################################################################
#### read files issued from the C++ phop algorithm #####
#if 'data' in filename : 

##    filename= 'data002'
#    data = np.loadtxt(filename)
#    atoms = data[:,0]
##    unique_times = data[:,1]
#    phops = data[:,2]
#    
#    plt.figure(1,[10,6])
#    base=1.1
#    Log_bins = np.array([p_threshold_low*base**i for i in range(-40, int(np.log(np.max(phops)/p_threshold_low)/np.log(base)), 1) ])
#    
##    for Tf in Tflist:
##        mask = ()
##        phops_tobin = phops[mask]
#    Nitem = len(phops)
#    P_phop, trash = np.histogram(phops, bins=Log_bins)
#    P_phop = P_phop / np.diff(Log_bins) / Nitem
#    plt.loglog(Log_bins[:-1], P_phop, lw=3, marker='o', label= '$p_{th}$='+filename[4:])
#    plt.legend(loc= 'best')
#    plt.xlabel(r'$p_{hop}$')
#    plt.ylabel(r'$N(p_{hop})$')
#    outName=filename + "_distro-phop"
#    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())





