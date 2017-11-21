## this computes:
#
import gsd.pygsd
import gsd.hoomd
import gsd.fl
import sys
from subprocess import call

outname = sys.argv[1] ## output name  (overwritten over)

MASTER_FILE='sources.py'
if len(sys.argv)>2: 
    MASTER_FILE = sys.argv[2]
files=[]
Master_Flow = open(MASTER_FILE,  'r')
for line in Master_Flow:
    if line[0] != '#':
        files.append(line[:-1])
Master_Flow.close()

print "merging all the following fiels into the output called '"+outname
for k in range(len(files)):
    print files[k]
print()

## copy #1 into #3
returnValue = call("cp "+files[0]+"  "+outname, shell=True)
#hoomdInputFlow1 = gsd.pygsd.GSDFile(open(filename1, 'rb'))
#hoomdTraj1 = gsd.hoomd.HOOMDTrajectory(hoomdInputFlow1)

## fail because not implemented
#hoomdOutputFlow3 = gsd.pygsd.GSDFile(open(outname, 'rw'))
#hoomdTraj3 = gsd.hoomd.HOOMDTrajectory(hoomdOutputFlow3)
#hoomdTraj3.extend(hoomdTraj2)

## open #2
for k in range(1,len(files)):
    with open(files[k], 'rb') as flow: 
        hoomdInputFlow2 = gsd.pygsd.GSDFile(flow)
        hoomdTraj2 = gsd.hoomd.HOOMDTrajectory(hoomdInputFlow2)
        print 'there are ', len(hoomdTraj2) , ' frames to be added  in file', files[k]
            
        ## open #3 in writing, append trajectory#2 to it.
        with gsd.fl.GSDFile(outname, 'wb') as f:
            t = gsd.hoomd.HOOMDTrajectory(f);
            t.extend( (hoomdTraj2[i] for i in range(len(hoomdTraj2)) ) )    ## creates an iteator, i.e. not a full copy of the flow !
        hoomdInputFlow2.close()

#print 'Done: file ', files[0] , '+', files[1], '  merged into  ', outname

