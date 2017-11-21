## this computes:
#   F_k(t), the so-called self-intermediate scattering function
#   MSD, 

import gsd.pygsd
import gsd.hoomd
import numpy as np
import sys
import module_filenamesHandler 
import module_computeObservables
import module_PBC

######
## plots: 
import module_importPlotParams
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True          
plt.rcParams.update(module_importPlotParams.fparams())
######
plt.ion()   ## optional

#TODO : check that selectedAtoms has been correctly used evreywhere

## DEBUG: with ovito  (hack a bit the cut of the files or reduce the diameters
#with open(kind+'_Ovito-dump.xyz', 'w') as flow:
#    flow.write(str(Natoms)+'\n')
#    flow.write("comment\n")
##    xyzAtom='C'
#    np.savetxt(flow, np.insert(positions,0,atomTypes,1), fmt='%d  %10.5f %10.5f %10.5f ')
##

def k_vector_norm(NX,NY,NZ, Lx):
    return np.sum( np.abs( 2j * np.pi/Lx * np.array([NX,NY,NZ]) )**2 )**0.5

################################################################################
####### Parameters selection ################
#############################################
## about time scales:
print 'usage:'
print 'run  compute_OBS*.py  [trajectory file.gsd]  [every_forMemory chosen=1]  [kind of the particles=A]    [dr=0.01]    [Nwindow_max=100]   [window_duration=100]'
filename = sys.argv[1]

every_forMemory = 1
if len(sys.argv)>2:
    every_forMemory = int(sys.argv[2])
# using  every_forMemory>1  is really just like increasing 
# recPeriod by a factor every_forMemory, it will just ignore many frames

kind = 'bulks8'     ## takes the bulks of the 8 cubes 
kind = 'bounds8'    # takes the boudaries of an 8-cube
kind = 'all'
kind = 'B'
kind = 'A'
if len(sys.argv)>3:
    kind = str(sys.argv[3])

dr = 0.01
if len(sys.argv)>4:
    dr = float(sys.argv[4])

Nwindow_max = 100
if len(sys.argv)>5:
    Nwindow_max = int(sys.argv[5])
    
window_duration = 100
if len(sys.argv)>6:
    window_duration = int(sys.argv[6])

################################################################################
####### read the file (using the gsd pakage) ##################################
with open(filename, 'rb') as flow:
    HoomdFlow = gsd.pygsd.GSDFile(flow)
    hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
    boxParams=(hoomdTraj.read_frame(0)).configuration.box
    Lx=boxParams[0]
    if boxParams[0] != boxParams[1] or  boxParams[1]!= boxParams[2] or  boxParams[0]!= boxParams[2]:
        print 'box dimensions are : ', boxParams[0]
        print 'and are not supported (only isotropic systems supported)'
        raise SystemExit
    size_factor = (int(Lx/9.3))
#    size_factor = (int(Lx/(9.3/2)))    ## debug only 
    Nframes = len(hoomdTraj)
    print 'there are Nframes=', Nframes, 'in the file, but we only use trajDuration=',Nframes/every_forMemory, ' of them.'
    Natoms = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'N')))
    expectedDensity = 1.2
    selectedAtoms, Nselected, expectedDensity = module_filenamesHandler.selectAtoms_function(hoomdTraj[0].particles.typeid, kind, Lx, expectedDensity, hoomdTraj[0].particles.position)
    trajDuration = Nframes/every_forMemory
    trajectory = np.zeros((trajDuration,Natoms,3))
    for time in range(0, trajDuration, 1):  
        ## we only look 1 frame every "every_forMemory" frame, to save memory: 
        trajectory[time] = (hoomdTraj[time*every_forMemory].particles.position) # [selectedAtoms]
    HoomdFlow.close()
    
print 'shape of your trajectory array:', np.shape(trajectory)
## "trajectory" is now the full array of all positions of all atoms (over time)
## trajectory[:,1] is the trajectory of particle 1 over time [its  (x,y,z) coords t be precise]
################################################################################


################################################################################
####################
### parameters: ####
####################

compute_RMSD=True
compute_RMSD=False

compute_Fkt=False
compute_Fkt=True

compute_gofr_Sq=True
compute_gofr_Sq=False

rMax = 9.5 # Lx/2. # 2.5 # Lx/2.0 ## -1 ## autoselect
shift_vector = np.array([1,0,0])*Lx*1.0/size_factor # g(r) shifted, to detect identical copies
#shift_vector = np.array([0,0,0])*0.    ## normal g(r) function
dt=0.0025

####### (DELTA) TIME RANGES USED #########
Delta_t_range = np.array([10,100])
Delta_t_range = np.array([1,3,10,30,100,300,1000,3000,10000,30000,100000])
Delta_t_range = np.array([0,1,2,3,6,10,30,60,100,300,600,1000,2400])
Delta_t_range = np.array([0,1,2,3,6,10,30,60,100,300,600,1000,3000, 6000, 10000, 30000, 100000, 300000,600000, 1000000])
#Delta_t_MAX = np.max(Delta_t_range)    # taking the max means all windows include all Delta_t's of the range.
Delta_t_MAX = 100     ## the Delta_t's larger than this value will not be computed for the latest windows.
## The last window will compute all Delta_t up to that one, Delta_t_MAX.
## one should have very long trajectories and then set this Delta_t to np.max(Delta_t_range)

## k vector used for F_k,s(t): #
NNN_sets = [ (4,6,8)] 
#NNN_sets = [ (1,0,0)] 

#############################################
####### End of parameters selection #########
################################################################################

recPeriod = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'rP')))      ## auto-detects the value of "recPeriod"
##############################
### simulation parameter (must be the same as what produced trajectory.gsd"
print 'Assuming that the dt in the MD ingration was: ', dt
##############################
#steps_per_tau_Hoomd =  1./dt
#steps_per_tau_KA    = (1./dt)/ 48**0.5
tau_KA_per_frame = recPeriod * dt # * 48**0.5       ## XXX : include this sqrt(48) iff you want to use the old standard of Kob and Andersen.
rootname_OBS=filename[:-4]+"_OBS="

Ws, every_forCPU, window_duration   =    module_computeObservables.makeWs_function(Delta_t_range, trajDuration, Delta_t_MAX, window_duration, Nwindow_max)

print('every_forCPU, window_duration, Ws = ', every_forCPU, window_duration, Ws)

################################################################################
### Fk(t) parameters (can be chosen reasonably) #####
#if compute_Fkt==True: 

print("choice of the direction of the vector k: ")
#NNN_sets = [ (k,k+2,k+6) for k in range(0,20,5)]
#NNN_sets = [ (10,10,16), (1,1,2),(1,1,5),(1,1,7), (2,2,5), (2,2,9), (2,2,16), (4,4,10),(4,4,1),(4,4,5),(4,4,12),(6,6,6),(6,6,9), (6,6,15), (8,8,10), (10,10,16),(10,10,26),(20,20,20)]    #NX -=1   :  so that we sample more differnet k's 
(NX,NY,NZ) = NNN_sets[0]

## Note that NNN_sets is then multiplied by this:
#    size_factor = 1
## so that larger systems are easily 
## studied using the approporiate equivalent set of k_vectors    
NX *= size_factor
NY *= size_factor
NZ *= size_factor
print '\nWe now use k=(', NX, NY, NZ, ')\n'
k_norm = k_vector_norm(NX, NY, NZ, Lx)
print '\nnorm of k vector:' , k_norm , '   (NX,NY,NZ)=', (NX,NY,NZ)
Fk_Dt      = np.zeros( len(Delta_t_range), dtype=complex)
Fk_Dt_std  = np.zeros( len(Delta_t_range), dtype=complex)
Fk_Deltat  = np.zeros( Natoms            , dtype=complex)
################################################################################

#if compute_RMSD==True:    
RMSD_t     = np.zeros( len(Delta_t_range) )
RMSD_t_std = np.zeros( len(Delta_t_range) )
RMSD_Deltat= np.zeros( Natoms )

if compute_Fkt or compute_RMSD :

    ## we compute things in windows of the time t0: [0, X], then [X,2X], etc. ##
    for win_index in range(len(Ws)-1):
        windowTag = "_wDur="+str(window_duration)+"_wBeg="+str(Ws[win_index]*recPeriod*dt)
        print '\nWe now compute window #', win_index

        ## we compute things for various time intervals $\Delta t$  ##
        Delta_t_index  = -1
        for Delta_t in Delta_t_range:
            Delta_t_index += 1

            Delta_t_weight = 0
            
            if compute_Fkt==True: 
                Fk_Deltat  *= 0.0

            if compute_RMSD==True:    
                RMSD_Deltat *= 0.0
            
            ## we compute the time-average of things, which at equilibrium is 
            ## equivalent to ensemble averaging.
            ## We do it "in parallel" for all particles, then average between particles
            interval_beg = Ws[win_index]
            interval_end = Ws[win_index+1]
            if (interval_end-Ws[win_index])/every_forCPU > 200 : 
                interval_end = interval_beg + every_forCPU
            for t0 in range(interval_beg,interval_end, every_forCPU):
                if t0+Delta_t > trajDuration-1 :
                    print "reached the end of this window's possible t_0's", t0, Ws[win_index]
                    break
                Delta_t_weight += 1
#                print("interval_beg,interval_end, every_forCPU",interval_beg,interval_end, every_forCPU)

                ## "displacements" includes all particles selected (e.g. type 'A')
#                [selectedAtoms]    TODO : adjsut using [selectedAtoms]
                displacements = np.abs(trajectory[t0+Delta_t,:]-trajectory[t0,:])     
                
                for t1 in range(len(displacements)):
                    displacements[t1] = module_PBC.applyPBC_vector_func(displacements[t1], Lx)
                
                if compute_Fkt==True:
                    ################  CORE  ###################
                    Fk_Deltat += module_computeObservables.Fkt_function(NX, NY, NZ, Lx, displacements)
                    ################  CORE  ###################

                if compute_RMSD==True:    
                    tmp = displacements**2
                    RMSD_Deltat += (tmp[:,0]+tmp[:,1]+tmp[:,2])**0.5    ## XXX is there really a sqrt here ??


            if Delta_t_weight > 0:
                if compute_Fkt==True:
                    Fk_Deltat /= (Delta_t_weight*6*2*2*2.0)
                    Fk_Dt[    Delta_t_index] = np.mean(Fk_Deltat)
                    Fk_Dt_std[Delta_t_index] = np.std( Fk_Deltat)  # variance wrt particles
                    print Delta_t*tau_KA_per_frame*every_forMemory, np.real(Fk_Dt[Delta_t_index]), np.real(Fk_Dt_std[Delta_t_index]), Delta_t
                if compute_RMSD==True:    
                    RMSD_Deltat /= Delta_t_weight
                    RMSD_t[Delta_t_index]    = np.mean(RMSD_Deltat)
                    RMSD_t_std[Delta_t_index]= np.std(RMSD_Deltat)  # variance wrt particles
            else:
                print 'The Delta_t=',Delta_t,' was not computed properly, the window was empty (or incomplete?)'
                break
        ### end of the loop on t0 ###
        #############################
            
        ## saving the F_k_t file:
        if compute_Fkt==True:
            tobesaved = np.zeros((len(Delta_t_range), 4))
            tobesaved[:,0] = Delta_t_range*tau_KA_per_frame*every_forMemory # deltat*recPeriod* dt * 48**0.5*every_forMemory = steps*dt * sqrt(48)
            tobesaved[:,1] = np.real(Fk_Dt)
            tobesaved[:,2] = np.real(Fk_Dt_std)
            tobesaved[:,3] = np.array(Delta_t_range , dtype=int)
            nonzero_indices = Fk_Dt!=0
            tobesaved = tobesaved[nonzero_indices]

            k_norm=k_vector_norm(NX, NY, NZ, Lx)
            k_norm_str = str(np.round(k_norm,2))
            outName=rootname_OBS+"Fkt_prt="+kind+"_k="+k_norm_str+"_triplet="+str(NX)+"_"+str(NY)+"_"+str(NZ)+'_samp='+str(Delta_t_weight)+windowTag
            
            with open(outName, 'a') as flow:
                header="$F_k(Delta t)$ for a bunch of  Delta_t values, using k=("+str(NX)+","+str(NY)+","+str(NZ)+")*2 pi/Lx , i.e. |k|="+k_norm_str+"    and for a window of time: "+ windowTag + "_wEnd(LJ units)="+str(Ws[win_index+1]*recPeriod*dt)+" , assuming dt="+str(dt)+"\n"
                header+= "# dt0 == "+str(every_forCPU*every_forMemory)+" (frames) == " +str(every_forCPU*every_forMemory*dt)+" (MD time units) == "+str(every_forCPU*every_forMemory*48**0.5)+" (KA time units) \n"
                header+= "Columns:  [ Delta time = Dt = Delta_Steps*sqrt(48)*dt , assuming dt="+str(dt)
                header+=  " ]   [ <F_k(Dt)> ]  [ stdDev(F_k(Dt)) ]   [ Delta_Steps ] "
                np.savetxt(flow, tobesaved, fmt="%f %f %f %d", header=header)
            
            ######### 
            tobesaved = np.loadtxt(outName)
            plt.figure(31,[10,6])
            plt.semilogx(tobesaved[:,0],tobesaved[:,1])
            ## TODO: add the error bars here, using stdDev to have their size.
            plt.title(r'Self-intermediate scattering')  #(self-spatial, 2-point temporal correlation function)
            plt.xlabel('$\Delta t$')
            plt.ylabel(r'$F_k(\Delta t)$')
            plt.ylim([0,1])
            plt.xlim([0.01,1e9])
            module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
            plt.close(31)

            
        ## saving the RMSD file: ##
        if compute_RMSD==True:    
            tobesaved = np.zeros((len(Delta_t_range), 4))
            tobesaved[:,0] = Delta_t_range*tau_KA_per_frame*every_forMemory
            tobesaved[:,1] = RMSD_t
            tobesaved[:,2] = RMSD_t_std
            tobesaved[:,3] = Delta_t_range
            nonzero_indices = RMSD_t!=0
            tobesaved = tobesaved[nonzero_indices]
        
            outName=rootname_OBS+"RMSD_prt="+kind+windowTag
            with open(outName, 'a') as flow:
                header ="$root Mean squared displacement (RMSD) for a bunch of  Delta_t values, for a window of time: "+ windowTag 
                header+= "_wEnd(LJ units)="+str(Ws[win_index+1]*recPeriod*dt)+" , assuming dt="+str(dt)+" \n"
                header+= "# dt0 == "+str(every_forCPU*every_forMemory)+" (frames) == " +str(every_forCPU*every_forMemory*dt)+" (MD time units) == "+str(every_forCPU*every_forMemory*48**0.5)+" (KA time units) \n"
                header+= "Columns:  [ Delta time = Dt = Delta_Steps*sqrt(48)*dt , assuming dt="+str(dt)
                header+= " ]   [ <RMSD> ]  [ stdDev(RMSD) ]   [ Delta_Steps ] "
                np.savetxt(flow, tobesaved, fmt="%f %f %f %d", header=header)

            tobesaved = np.loadtxt(outName)
            plt.figure(30,[10,6])
            plt.semilogx(tobesaved[:,0],tobesaved[:,1])
            ## TODO: add the error bars here, using stdDev to have their size.
            plt.title(r'Root Mean Square Displacement, RMSD$(\Delta t)$')
            plt.xlabel('$\Delta t$')
            plt.ylabel(r'RMSD$(\Delta t)$')
            plt.ylim([0,2.0])
            plt.xlim([0.01,1e9])
            module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
            plt.close(30)

## : ##
if compute_gofr_Sq==True:    
    ## we load a single Frame:
#    FrameToChoosePositionsFrom = len(trajectory)-40
#    print "We compute G(r) and S(q), using the frame # ", FrameToChoosePositionsFrom

    print ("Computing g(r) using a shift vector (to detect decorrelation between identical sub-samples within the sample)  of shift_vector=", shift_vector, "   and a size_factor (autodetect multiplicative factor of the initial system) of size_factor=", size_factor , " within the windows:" , Ws)

    ## we compute things in windows of the time t0: [0, X], then [X,2X], etc. ##
    base = 1.6
    win_index_log_max = int(np.log(len(Ws)  )/np.log(base))+2
    win_index_max = len(Ws)-1
    if Nwindow_max > 1000: 
        ## activate the log-scaling of the windows #
        win_index_max = win_index_log_max
    for kpow in range(0, win_index_max):
        if win_index_max == win_index_log_max :
            ## use the log-scaling of the windows #
            win_index = int(base**kpow)
        else:
            win_index = kpow
    
        print '\nWe now compute window #', win_index

        Sq_average = 0.0
        gofr_average  = 0.0
        Delta_t_weight = 0.0
        
#        every_forCPU=window_duration/1-1    ## XXX dirty hack to hae much less t0 's sampled 
#        for t0 in range(1): #Ws[win_index],Ws[win_index+1], every_forCPU):
        N_snapshot_to_average_over = min(5, window_duration)
        first_snapshot_to_average_over = Ws[win_index]

        last_snapshot_to_average_over = min(Ws[win_index+1], Ws[win_index]+N_snapshot_to_average_over)
        every_snapshot_to_average_over = 1

#        last_snapshot_to_average_over = Ws[win_index+1]
#        every_snapshot_to_average_over = max(1, (last_snapshot_to_average_over -  first_snapshot_to_average_over) / (N_snapshot_to_average_over+1))

        N_samples_averaged_over = 0
        last_snapshot_to_average_over = min(last_snapshot_to_average_over, trajDuration)
        for t0 in range( first_snapshot_to_average_over, last_snapshot_to_average_over, every_snapshot_to_average_over ): 
            N_samples_averaged_over += 1

            FrameToChoosePositionsFrom = t0 
            positions = trajectory[FrameToChoosePositionsFrom]
            rvalues, gofr_normalized, qs, Sq =  module_computeObservables.gofr_and_Sq_computer(positions, selectedAtoms, dr, Lx, expectedDensity, shift_vector, rMax)
            print("function call finished")
    
            gofr_average += gofr_normalized
            Sq_average += Sq
            Delta_t_weight +=1.0

        gofr_normalized = gofr_average / Delta_t_weight
        Sq = Sq_average / Delta_t_weight
        
        ## TODO: smooth g(r) using sliding windows, so that high-freq stuff go away? 
        ## it would be more refined than just taking larger dr, yet ... not really worth it
        
        tobesaved = np.zeros((len(gofr_normalized), 4))
        tobesaved[:,0] = rvalues
        tobesaved[:,1] = gofr_normalized
        tobesaved[:,2] = qs
        tobesaved[:len(Sq),3] = Sq
    #    nonzero_indices = gofr_normalized!=0
    #    tobesaved = tobesaved[nonzero_indices]

#        windowTag = "_wDur="+str(window_duration)+"_wBeg="+str(Ws[win_index])
        windowTag ='_samp='+str(N_samples_averaged_over)+"_wBeg="+str(Ws[win_index]*recPeriod*dt)   
        outName=rootname_OBS+"gofr_dr="+str(dr)+"_rMax="+str(np.round(rMax,1))+"_prt="+kind+ windowTag #"_Fr="+str(FrameToChoosePositionsFrom)
#        +

        with open(outName, 'w') as flow:
            header ="$ g(r) and its fft S(q) for the frame number "+str(FrameToChoosePositionsFrom)+"\n"
            header+= "Columns:  [ r ]   [ g(r) ]  [ q ]   [ S(q) ] "
            np.savetxt(flow, tobesaved, fmt="%f %f %f %f", header=header)

        tobesaved = np.loadtxt(outName)
        plt.figure(33,[10,6])
        plt.plot(tobesaved[:,0],tobesaved[:,1])
        plt.xlabel('$r$')
        plt.ylabel(r'$g(r)$')
        plt.ylim([0,2.5])
        module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
        plt.close(33)
        
        plt.figure(34,[10,6])
        plt.plot(tobesaved[1:,2],tobesaved[1:,3])
        plt.xlabel('$q$')
        plt.ylabel(r'$S(q)$')
        plt.xlim([-5,60])
    #    plt.ylim([0,4])
        module_importPlotParams.savefig_perso(outName+'Sq'+module_importPlotParams.fsaveFigFormat())
        plt.close(34)


#    L=len(gofr_normalized)
#    a=np.zeros(2*L)
#    a[:L] = gofr_normalized
#    a[L:] = gofr_normalized[::-1]
    

#        else: 
#            pass
#            plt.xlim([0.1,50])
##            tobesaved[NNN_sets_index,3] = k_norm
#            flow=open(rootname_OBS+"Fkt_prt="+kind+"_Deltat="+str(Delta_t_range[0]), 'a')  
#            module_importPlotParams.savefig_perso(rootname_OBS+"Fkt_prt="+kind+"_Deltat="+str(Delta_t_range[0])+module_importPlotParams.fsaveFigFormat())
#            flow.write()   ## 
###        np.savetxt(flow, tobesaved)
#        flow.close()
    #        plt.close(31)




