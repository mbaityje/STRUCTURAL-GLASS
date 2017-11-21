## B1-construct_TrainingSet.py
## Construct the training set for the Machine Learning Method
##    Essentially just pick a few examples from the windows of very still and 
##    very fast particles time-ranges that were detected in A2.py 
##
## we then actually go and read the file of (FIREd) trajectory, called trajname. We go only to the "few" 2*Nexamples atoms that we need for training, and compute+record their descriptors to a file.
## 
## input:   trainers-hard.dat  
##          trainers-soft.dat
##          FIRE-trajectory*.gsd  OR  decorated-trajectory*.gsd   (to be coded, though.. taking care of reading the correct time steps)
## output: screen+
##          hard-soft.gz ?

from __future__ import print_function
import gsd.pygsd
import gsd.hoomd
import gsd.fl

import sys
import numpy as np

import module_filenamesHandler 
import module_PBC
import module_descriptors

######
## plots: 
import module_importPlotParams
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True          
plt.rcParams.update(module_importPlotParams.fparams())
plt.ion()   ## optional
######

'''
Notes: 
One may define hard particles in many ways:
Ask for a slowWindow=500  (size)
- take the half
- take the state after some time of inactivity (maximize the time before an event)
- take 100 steps before next activation 
- +take particles that never moved, ever.. (in our sample)

Soft particles:
- find slow intervals, and take a few steps before
- use only the "epicenters" of avalanches

'''

print('usage:')
print('run  '+sys.argv[0]+'   [rootname_AVA]   [Nexamples]   ')

verbose = True
verbose = False

rootname_AVA =     sys.argv[1]
Nexamples     = int(sys.argv[2])

throwBoundaries_steps = 10
soft_delay = 2 ## time frames, which is not a well controlled unit

### parameters to define hard and soft instants ###
#hard_factor = 2 

Natoms = int(float(module_filenamesHandler.filename_parser(rootname_AVA[:-4], 'N')))
p_threshold_low = max(float(module_filenamesHandler.filename_parser(rootname_AVA[:-4], 'pl')), float(module_filenamesHandler.filename_parser(rootname_AVA[:-4], 'pLow')))
p_threshold_high=float(module_filenamesHandler.filename_parser(rootname_AVA[:-4], 'ph'))
slowWindow = int(float(module_filenamesHandler.filename_parser(rootname_AVA[:-4], 'sW')))
fastWindow = int(float(module_filenamesHandler.filename_parser(rootname_AVA[:-4], 'fW')))
kind                 = module_filenamesHandler.filename_parser(rootname_AVA[:-4], 'kind')
curStep   = int(float(module_filenamesHandler.filename_parser(rootname_AVA[:-4], 'cst'))/1e6)
recPeriod = int(float(module_filenamesHandler.filename_parser(rootname_AVA[:-4], 'rP')))
LL = int(float(module_filenamesHandler.filename_parser(rootname_AVA[:-4], 'LL')))

rootname = module_filenamesHandler.get_rootname(rootname_AVA)
suffix   = module_filenamesHandler.get_suffix(rootname_AVA)
softwindowsName = rootname_AVA + "trainers-soft.dat"
hardwindowsName = rootname_AVA + "trainers-hard.dat"
trajname = rootname+"_type=FIRE-traj_"+"cst="+str(curStep)+"e6_rP="+str(recPeriod)+"_LL="+str(LL)+".gsd"
typeTraj = module_filenamesHandler.filename_parser(trajname[:-4], 'type')
if "FIRE" not in typeTraj or 'deco' in typeTraj: 
    print("error: use FIRE trajectories (in their eventsData form) OR adapt your code")
    print("\nNOTE: we need the full trajectory when we fish for the ebst exampes (to train our SVM), because we will use some frames that are in the middle of nowhere, and also just a few steps BEFORE events. So the 'deco'  file is not enough.")
    raise SystemExit



print("you selected only atoms of the kind ", kind)
################################################################################
#### reading the "decorated" file, and plot stuff #####
if "AVA=" in rootname_AVA :   

#note: le but de ce code est maintenant juste de lire le nombre d'intervalles med, soft et hard, 
#et de tirer de chaque intervalle hard/soft 1 et 1 seule instant d'interet. 
#il ne ffaut pas faire de re-selection ici, on doit se contenter de lire les intrervalles deja tries sur le volet, et de piocher 1 config par intervalle.

    softwindows = np.loadtxt(softwindowsName, dtype=int)
    hardwindows = np.loadtxt(hardwindowsName, dtype=int)
    

#    outname = eventsDataname[:-4]+"_slowW="+str(slowWindow)+"_part=hard-soft"
#    np.savetxt(outname, AtomTimesAll[:,:2], fmt="%d %d", header="hard, THEN soft particles found using slowWindow as in title. \nColumns: atom number    time step ") 
#    print("The training instants have been shuffled and recorded to the file:\n",outname,'\n')
    
    
    ### fish out a few instants from rthe windows ####
    print("We now fish (fetch) out one instant per window.")
    AtomTimesSoft = np.zeros((len(softwindows),2),dtype=int)
    soft_num=0
    for i in range(len(softwindows)):
        beg = softwindows[i,1]
        if beg-soft_delay > throwBoundaries_steps :
            AtomTimesSoft[soft_num,0]= softwindows[i,0]
            AtomTimesSoft[soft_num,1]= beg - soft_delay     ## taking a few steps before the begginging of the fastWindow
            soft_num+=1
    AtomTimesSoft = AtomTimesSoft[:soft_num]

    AtomTimesHard = np.zeros(((len(hardwindows),2)),dtype=int)
    hard_num=0
    for i in range(len(hardwindows)):
        beg = hardwindows[i,1]
        end = hardwindows[i,2]
        if beg > throwBoundaries_steps and end < LL - throwBoundaries_steps :
            AtomTimesHard[hard_num,0]= hardwindows[i,0]
            AtomTimesHard[hard_num,1]= (beg+end)/2      ## taking the middle of the slowWindow
            hard_num+=1
    AtomTimesHard = AtomTimesHard[:hard_num]

    n = min(soft_num, hard_num)
    if n < Nexamples :
        Nexamples = n
        print("You don't have that much data. I restrict you to only Nexamples=", Nexamples)
    del n
    
    #####
    ### shuffle the instants before cropping them ###
    seed_shuffle = 42
    np.random.seed(seed_shuffle)

    soft_shuffle = np.arange(Nexamples, dtype=int)
    np.random.shuffle(soft_shuffle) ## in-place operation
    AtomTimesSoft = AtomTimesSoft[soft_shuffle[0:Nexamples],:]

    hard_shuffle = np.arange(Nexamples, dtype=int)
    np.random.shuffle(hard_shuffle) ## in-place operation
    AtomTimesHard = AtomTimesHard[hard_shuffle[0:Nexamples],:]

    AtomTimesAll          = np.zeros((2*Nexamples, 2), dtype=int)
    AtomTimesAll[ 0      :  Nexamples] = AtomTimesHard
    AtomTimesAll[Nexamples:2*Nexamples] = AtomTimesSoft



    print("Data has been fished, shuffled and cut. Now reading the trajectory file to get positions.")    
    filename = trajname
    every_forMemory=1
    unique_times = np.sort(np.array(list( set( AtomTimesAll[:,1]) )))  #+ list(AtomTimesHard[:,1])
    descriptedDuration = len(unique_times)
    print("There is a number ", descriptedDuration, " of unique times(=frames) to look at")
    backArray = np.zeros(np.max(unique_times)+1, dtype=int)-1
    for tc_index in range(descriptedDuration):
        backArray[unique_times[tc_index]] = tc_index
    ## we now have:    backArray[unique_times[:]]  is range(descriptedDuration)
    ## and backArray[AtomTimesAll[:,1]] is the stuff we'll use
    
    ################################################################################
    ####### read the file (using the gsd pakage) ##################################
    with open(filename, 'rb') as flow:
        HoomdFlow = gsd.pygsd.GSDFile(flow)
        hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
        atomTypes=np.array((hoomdTraj.read_frame(0)).particles.typeid, dtype=int)
        boxParams=(hoomdTraj.read_frame(0)).configuration.box
        Lx=boxParams[0]
        if boxParams[0] != boxParams[1] or  boxParams[1]!= boxParams[2] or  boxParams[0]!= boxParams[2]:
            print('box dimensions are : ', boxParams[0])
            print('and are not supported (only isotropic systems supported)')
            raise SystemExit
        Nframes = len(hoomdTraj)
        
#        expectedDensity = float(module_filenamesHandler.filename_parser(filename[:-4], 'rho'))
#        selectedAtoms, Natoms, expectedDensity = module_filenamesHandler.selectAtoms_function(atomTypes, kind, Lx, expectedDensity, 0)
        
        trajectory = np.zeros((descriptedDuration, Natoms, 3))
        print('there are Nframes=', Nframes, 'in the file, but we only use trajDuration=',descriptedDuration, ' of them.')
        for tc_index in range(descriptedDuration):  
            interesting_time = unique_times[tc_index]
            ## we only look 1 frame every "every_forMemory" frame, to save memory: 
            trajectory[tc_index] = (hoomdTraj[interesting_time].particles.position)  ## we DO NOT use selectedAtoms because we need all atoms to compute the descriptor functions
        HoomdFlow.close()
        
    print('shape of your trajectory array (interesting_time frames only):    ', np.shape(trajectory))
#    print('shape of your trajectory array (non empty): ', np.shape(trajectory[unique_times]))
    







    ###############################################
    ##### THE CORE COMPUTATION HAPPENS HERE : #####
    ## computing the structure functions ##
    print("\nStarting the conversion of instants  (atom, time)  into descriptors  (atom, time, {G_i}_i)")
    descriptor_family_parameters = np.loadtxt("B0_fam=1_type=radial.dat")
    descriptor_family_name=1
    all_descriptors_array = module_descriptors.descriptor_func_allInstants(AtomTimesAll, trajectory, backArray, descriptor_family_name, descriptor_family_parameters, Lx, atomTypes)
    ###############################################
    
    
    
    
    
    

    outname = rootname_AVA+"descriptors_fam="+str(descriptor_family_name)+"_Nex="+str(Nexamples)+".dat"
    header = "hard, THEN soft particles found using slowWindow as in title. Each line corresponds to an instant. Each column is a {G_i} value. "
    np.savetxt(outname, all_descriptors_array, header=header) 
    print("descriptors list saved to file(s) (compressed and not compressed):\n", outname, "\n" )

#    ## compressed version (can still be opened with   np.loadtxt(name)  ) ##
#    outname = ... +".gz"
#    with open(outname, 'wb') as flow:
#        np.savez(flow, all_descriptors_array, header=header) 




