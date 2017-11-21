## this computes:
#
import gsd.pygsd
import gsd.hoomd
import gsd.fl
import sys
from subprocess import call

filename1 = sys.argv[1]
filename2 = sys.argv[2]
filename3 = sys.argv[3] ## output name  (overwritten over)


## copy #1 into #3
returnValue = call("cp "+filename1+"  "+filename3, shell=True)
#hoomdInputFlow1 = gsd.pygsd.GSDFile(open(filename1, 'rb'))
#hoomdTraj1 = gsd.hoomd.HOOMDTrajectory(hoomdInputFlow1)

## fail because not implemented
#hoomdOutputFlow3 = gsd.pygsd.GSDFile(open(filename3, 'rw'))
#hoomdTraj3 = gsd.hoomd.HOOMDTrajectory(hoomdOutputFlow3)
#hoomdTraj3.extend(hoomdTraj2)

## open #2
with open(filename2, 'rb') as flow: 
    hoomdInputFlow2 = gsd.pygsd.GSDFile(flow)
    hoomdTraj2 = gsd.hoomd.HOOMDTrajectory(hoomdInputFlow2)
    print 'there are ', len(hoomdTraj2) , ' frames to be added '    
        
    ## open #3 in writing, append trajectory#2 to it.
    with gsd.fl.GSDFile(filename3, 'wb') as f:
        t = gsd.hoomd.HOOMDTrajectory(f);
        t.extend( (hoomdTraj2[i] for i in range(len(hoomdTraj2)) ) )    ## creates an iteator, i.e. not a full copy of the flow !
        
    print 'Done: file ', filename1 , '+', filename2, '  merged into  ', filename3

    hoomdInputFlow2.close()

