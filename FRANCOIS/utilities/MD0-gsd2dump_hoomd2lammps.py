## this computes:
#
import gsd.pygsd
import gsd.hoomd
import numpy as np
import sys
from subprocess import call 
import module_filenamesHandler


filename = sys.argv[1]
every=1


folder= "States-"+filename[:-4]+ '/'
returnValue = call("mkdir "+folder , shell=True)  ## !! option shell=True is quite dangerous !!



############################
### read the file (gsd) ####
with open(filename, 'rb') as flow : 
    hoomdInputFlow = gsd.pygsd.GSDFile(flow)
    hoomdTraj = gsd.hoomd.HOOMDTrajectory(hoomdInputFlow);
    print 'using hoomdTraj.file.gsd_version=', hoomdTraj.file.gsd_version
    print 'using hoomdTraj.file.gsd_application=', hoomdTraj.file.application
    print 'using hoomdTraj.file.schema=', hoomdTraj.file.schema
    Nframes = hoomdTraj.file.nframes

    ## we assume these do not change over time:
    boxParams   = hoomdTraj.read_frame(0).configuration.box
    dim         = hoomdTraj.read_frame(0).configuration.dimensions
    Natoms      = hoomdTraj[0].particles.N
    atomIds     = np.arange(Natoms, dtype=int).transpose()
    atomTypes   = np.array(hoomdTraj[0].particles.typeid, dtype=int).transpose()

    # precision can be set with this %.10f thingy
#    precision = "%d %d %.10f %.10f %.10f"
    precision = "%d %.10f %.10f %.10f"
    Lx = boxParams[0]
    if boxParams[0]==boxParams[1] and boxParams[1]==boxParams[2] :
        print 'ok, system is isotropic.'

        for t0 in range(0,Nframes,every):
            step = hoomdTraj.read_frame(t0).configuration.step
            if t0==0 or t0>Nframes-every:
                print "step transformed: ", step
            ## assuming boxParams does not change
            header = module_filenamesHandler.LAMMPS_header_function(step, Natoms, boxParams)
            
            pos = (hoomdTraj[t0].particles.position + Lx/2.0)/Lx
#            pos = np.insert(pos,  obj=0, values=atomTypes+1, axis=1)
            pos = np.insert(pos,  obj=0, values=atomIds+1, axis=1)
            with open(folder+"States0."+str(step)+'.dump', 'w') as outFlow:
                np.savetxt(outFlow, pos, fmt=precision, header = header, comments='')   
            
        hoomdInputFlow.close()

    else: 
        print 'non cubic box not supported (easy to add though)'
    
print 'Convert complete. Be careful that my atomTypes are indexed starting from 0, not from 1 !!'
    
    
    
