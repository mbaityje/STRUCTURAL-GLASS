# module_filenamesHandler.py
#import None

#rootname= "KA_rho=12e-1_N="+str(Natoms)+"_T="+str(temperatureTimes100)+"e-2_tag="+special_tag+"_type="
def log_action(rootname, action):
    print '[action logged:]\n'+action
    with open(rootname+'_type=notes.txt', 'a') as flow:
        flow.write(action)
    return 

def filename_parser(filename, variable):
    # for example to find the value of p1 in filename, put variable=='p1'
    variableCopy=variable+'='
    filename_splitted = filename.split('_')
    for word in filename_splitted:
        if variableCopy in word :
            return word[len(variableCopy):]
    return 'error: the keyword   '+variable+'   is not present'

## copy all words into rootname, until the word "type=" is reached
def get_rootname(filename):
    variable='type'
    variableCopy=variable+'='
    filename_splitted = filename.split('_')
    rootname=''
    for word in filename_splitted:
        if variableCopy not in word :
            rootname+=word+'_'
        else: 
            break
    return rootname[:-1]
    
    
## copy all words (from the end) into rootname, before the word "type=" is reached
def get_suffix(filename):
    variable='type'
    variableCopy=variable+'='
    filename_splitted = filename.split('_')
    suffix=''
    for word in filename_splitted[::-1]:
        if variableCopy not in word :
            suffix = word+'_'+suffix
        else: 
            break
    return suffix[:-1]

## copy all words into rootname, until the word "N=" is reached
def get_root_rootname(filename):
    variable='N'
    variableCopy=variable+'='
    filename_splitted = filename.split('_')
    rootname=''
    for word in filename_splitted:
        if variableCopy not in word :
            rootname+=word+'_'
        else: 
            break
    return rootname[:-1]


def LAMMPS_header_function(step, Natoms, boxParams):
    header  = "ITEM: TIMESTEP\n"
    header += str(step)+"\n"
    header+= "ITEM: NUMBER OF ATOMS\n"
    header += str(Natoms)+"\n"
    header+= "ITEM: BOX BOUNDS pp pp pp\n"
    header += str(0)+ " "+str(boxParams[0])+"\n"
    header += str(0)+ " "+str(boxParams[1])+"\n"
    header += str(0)+ " "+str(boxParams[2])+"\n"
    header+= "ITEM: ATOMS id xs ys zs"  ## this string needs to be EXACTLY one of the three types:
    ####string dumpfieldsscaled("ITEM: ATOMS id xs ys zs"); // style 0
    ####string dumpfieldscna("ITEM: ATOMS id x y z StructureType"); // style 1
    ####string dumpfieldsall("ITEM: ATOMS id xs ys zs xu yu zu c_csym c_free_energy c_sigma[1] c_sigma[2] c_sigma[3] c_sigma[4] c_sigma[5] c_sigma[6] c_c
    return header    

    
## this is a very light function, which selects populations of
## particles, in order to compute stuff only for this sub-set of particles
def selectAtoms_function(atomTypes, kind, Lx, expectedDensity, positions_ToChooseAtomsFrom):

    ## extract atomtypes from frame 0, since types do not change:
#    atomTypes = hoomdTraj[0].particles.typeid  
    expectedDensity /= len(atomTypes)
    validList = ['all', 'A', 'B', 'bulks8', 'bounds8' ]
    if kind not in validList:
        print "\nError: incorrect 'kind' argument. Use an existing sub-set of particles, i.e. one of these:  "+str(validList)+"\n"
        raise SystemExit
    
    if kind == 'all':
        selectedAtoms = (atomTypes != -1)   ## this is always true
    
    elif kind == 'A':
        selectedAtoms = (atomTypes == 0)
        
    elif kind == 'B':
        selectedAtoms = (atomTypes == 1)
        
    elif kind == 'bulks8' or kind == 'bounds8':
        ## we track particles according to where they were
        ## at this chosen time (FrameToChooseAtomsFrom) : 
#        positions = hoomdTraj[FrameToChooseAtomsFrom].particles.position
        positions = positions_ToChooseAtomsFrom
        
        ## divide space in bulk pieces and boudnary pieces, assuming 
        ## that it is made of 8 cubes of size L=Lx/2
        ## Each cube is divided in a central cube of size a^3=L^3/2
        ## plus its boudnaries of size L^3-a^3 = L^3/2
        ## thus a = L/2**0.5
        ## thus pieces with  -L+
        a = (Lx/2)/2**(1/3.)
        selectedAtoms = np.ones(len(positions), dtype=bool)
        for dim in range(3):
            selectedAtoms *= ( ((positions[:,dim] > -Lx/4-a/2 )*( positions[:,dim] < -Lx/4+a/2 )) +             ((positions[:,dim] >  Lx/4-a/2 )*( positions[:,dim] <  Lx/4+a/2 )) )
        ### DEBUG:    
        #print np.shape(positions), Lx, a, Lx/2
        #print -Lx/4-a/2, -Lx/4+a/2, Lx/4-a/2, Lx/4+a/2

    ## invert the selection ##
    if kind == 'bounds8':   
        ## atoms on the boundaries are exactly atoms not in the bulk 
        selectedAtoms = np.array(1.0-selectedAtoms*1.0, dtype=bool)

    ## morally, we now have: ##
    atomTypes = atomTypes[selectedAtoms]
    expectedDensity *= len(atomTypes)
    Natoms = len(atomTypes)
    return selectedAtoms, Natoms, expectedDensity
    


##steps_per_tau_Hoomd =  1./dt
##steps_per_tau_KA    = (1./dt)/ 48**0.5
###steps_per_tau_alpha =((1./dt)/ 48**0.5) * 1e4   # 288675.13459481293 = 3e5 if dt=0.005
###steps_per_tau_alpha = 6e5    ## NVT time steps
###record_time_nve = 100 * steps_per_tau_alpha ## NVE time steps

