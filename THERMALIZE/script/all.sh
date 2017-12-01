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
CheckThermalizationAll.sh

#Find the trajectory of the IS through a bisection algorithm.
#Launches ChunkSect.sh
ChunkSectAll.sh 
