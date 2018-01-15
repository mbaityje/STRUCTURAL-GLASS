exit #Evitemos cagadas

#Logical succession of all the scripts, to perform all the tasks in
#the repository.

# Create the initial configuration that will be fused a ton of
# times. Only needs to be called once.
bash CreateInitialIS.sh

#Launch thermalization runs for all temperatures
bash ThermalizeAll.sh

#Check thermalization for all the systems. Calculates msd, Fkt and tau.
#Launches SelfIntermediateScatteringFunction.sh
bash CheckThermalizationAll.sh

#Find the trajectory of the IS through a bisection algorithm.
#Launches ChunkSect.sh
bash ChunkSectAll.sh 

#From an initial thermal configuration calculates the overlap and the overlap of the IS (as well as the energies) as a function of time.
bash OverlapTrajectoryIS.sh

#Given a trajectory, reruns the minimization at all EVERY time step starting from the frame $iframe, and measure Eis, msd_IS of the time step, and number of particles that moved in that time step.
bash MinimizeSegment.sh



