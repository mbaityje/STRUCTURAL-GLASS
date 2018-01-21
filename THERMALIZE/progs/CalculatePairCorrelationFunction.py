###########################################################
#                                                         #
# This utility counts the number of frames in a gsd file. #
#                                                         #
###########################################################


from __future__ import print_function #for compatibility with python3.5
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
import hoomd #hoomd is the MD package from the Glotzer group
import argparse
from lib import module_measurements as med
from matplotlib import pyplot as plt



#Start hoomd
hoomd.context.initialize()
hoomd.option.set_notice_level(0)
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd



#Read parameters
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename', #positional argument
                    nargs=1,
                    help='name of the .gsd trajectory we want to read'
)
parser.add_argument('--dr', #optional argument
                    nargs=1,
                    type=float,
                    required=False,
                    default=[0.03],
                    help='dr for the binning of the g(r)'
)
parser.add_argument('--rMax', #optional argument
                    nargs=1,
                    type=float,
                    required=False,
                    default=[-0.5],
                    help='Maximum radius for the g(r). rMax>0: absolute value. rMax<0: we set it as \
                    a fraction of L. For example, rMax=-0.5 results in rMax=L/2.'
)
parser.add_argument('-l','--label', #optional argument
                    nargs=1,
                    required=False,
                    default=[''],
                    help='label for distinguishing runs and continuations'
)




args = parser.parse_args(more_arguments)
filename=args.filename[0]
dr=args.dr[0]
rMax_temp=args.rMax[0]
label=str(args.label[0])
del parser



print("filename: ",filename)
print("dr = ",dr)
print("rMax_temp = ",rMax_temp)
print("label = ",label)




system = hoomd.init.read_gsd(filename=filename)
Natoms = len(system.particles)
assert(system.box.Lx == system.box.Ly == system.box.Lz)
L=system.box.Lx
rMax=rMax_temp if rMax_temp>=0 else -rMax_temp*L
assert(0< rMax <= 0.5*L)
assert(system.box.get_volume()==L*L*L)
rho=Natoms/(system.box.get_volume())
print("Natoms = ",Natoms)
print("rho = ",rho)
print("rMax = ",rMax)



################################################################
# 
# CALCULATE PAIR CORRELATION FUNCTION
#
################################################################
snapshot=system.take_snapshot(dtype='double')
positions=np.array(snapshot.particles.position) #Now 'positions' is a NatomsxDIM vector,
                                                #storing all the particles' positions
distances=med.CalculateRelativeDistances(positions,Natoms,L) #A list of the Natoms*(Natoms-1)/2 relative distances
gofr=med.CalculatePairCorrelationFunction(distances, Natoms, dr=dr, rMax=rMax, number_density=rho)
rvalues = np.arange(0,rMax,dr)
output_gofr=np.column_stack((rvalues, gofr))
np.savetxt('gofr'+label+'.txt',output_gofr,fmt='%g %.14g')


#I plot the function I found
plt.plot(rvalues, gofr)
plt.xlabel('r')
plt.ylabel('g(r)')
plt.title('Pair Correlation function of '+filename)
plt.grid(True)
plt.savefig("gofr.png")

