### B3-assign-softness.py
from __future__ import print_function

import gsd.pygsd
import gsd.hoomd
import gsd.fl
import numpy as np
import sys

import module_filenamesHandler 
import module_descriptors
import module_dist_compute

print('usage:')
print('run  '+sys.argv[0]+'   [deco-traj.gsd]   [hyperplane_name.dat??] ')

# run   B3-assign-softness.py  KA_rho=12e-1_N=1e3_T=43e-2_type=FIRE-traj-deco_cst=91e6_rP=1000_LL=60000_SLID_pLow=0.002_rec=100.gsd     KA_rho=12e-1_N=1e3_T=43e-2_type=FIRE-traj-deco_cst=91e6_rP=1000_LL=60000_SLID_pLow=0.002_rec=100_kind=A_pl=0.002_ph=0.08_sW=500_fW=2_AVA=descriptors_fam=1_Nex=10000_hyperplane.dat  

filename        = sys.argv[1]
hyperplane_name = sys.argv[2]
outputName_softness_R=filename[:-4]+"_SOFT=softness-R.gsd"

Natoms = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'N')))
kind   =           module_filenamesHandler.filename_parser(hyperplane_name[:-4], 'kind')


descriptor_family_parameters = np.loadtxt("B0_fam=1_type=radial.dat")
descriptor_family_name=1

#    hyperplane_name = filename[:-4]+"_hyperplane.dat"
hyperplane = np.loadtxt(hyperplane_name)
descriptor_averages = hyperplane[:,0]
descriptor_stdDevs  = hyperplane[:,1]
a_hyperplane        = hyperplane[:,2]
with open(hyperplane_name, 'r') as flow:
    flow.readline()
    b_intercept = float((flow.readline()).split()[-1])



if "deco" in filename :   # read the decorated file produced before # 

    p_threshold_low = float(module_filenamesHandler.filename_parser(filename[:-4], 'pLow'))
    expectedDensity = float(module_filenamesHandler.filename_parser(filename[:-4], 'rho'))
    print('reading that the file has a p_threshold_low=', p_threshold_low, ' ...')
    if p_threshold_low ==0 : 
        p_threshold_low= 1e-7
        print("correcting its threshold to :",p_threshold_low)

    Natoms = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'N')) )


#    # TODO write a function where you open the output "decorated" file in writing to add jsut the softness thing to it (?)
#    ff = gsd.fl.GSDFile(name=filename, mode='wb')
#    decoratedFlowRead = gsd.fl.GSDFile(name=filename, mode='rb')
#decoratedFlowRead.read_chunk(frame=tc_index, name='softness')
    with gsd.fl.GSDFile(name=filename, mode='rb') as decoratedFlowRead :    # TODO : use with as instead of close() --- everywhere !!
        atomTypes = np.array(decoratedFlowRead.read_chunk(frame=0, name='atomTypes'), dtype=int)
        atomIds   = decoratedFlowRead.read_chunk(frame=0, name='atomIds')
        boxParams   = decoratedFlowRead.read_chunk(frame=0, name='configuration.box')   
        Lx = boxParams[0]

        gsd.fl.create(name=outputName_softness_R,application="trajectory with decorations",\
                schema="[atomTypes] [atomIds] tc positions phop TAs TBs averageA Prev_StdDev Next_StdDev average_distAB softness",\
                schema_version=[1,0])
        with gsd.fl.GSDFile(name=outputName_softness_R, mode='wb') as ff: 
            # writing atom types only in the first frame : (and not closing this first frame) #
            ff.write_chunk(name='atomTypes' , data=np.array(atomTypes , dtype=int))
            ff.write_chunk(name='configuration.box' , data=np.array(boxParams , dtype=np.float32))
            ff.write_chunk(name='atomIds' , data=np.array(atomIds , dtype=int))


            unique_times_number = decoratedFlowRead.nframes
            unique_times    = np.zeros(unique_times_number, dtype=int)
    #        trajectory      = np.zeros((unique_times_number, Natoms, 3))
            phop           = np.zeros((unique_times_number, Natoms))
        #    times           = np.zeros((unique_times_number, Natoms), dtype=int)
        #    TAs             = np.zeros((unique_times_number, Natoms), dtype=int)
        #    TBs             = np.zeros((unique_times_number, Natoms), dtype=int)
        #    averageA        = np.zeros((unique_times_number, Natoms, 3))
        #    Prev_StdDev     = np.zeros((unique_times_number, Natoms))
        #    Next_StdDev     = np.zeros((unique_times_number, Natoms))
        #    average_distAB  = np.zeros((unique_times_number, Natoms))
            softness        = np.zeros((unique_times_number, Natoms))
            for tc_index in range(3): # unique_times_number): ## TODO : oput all tiems here 
                unique_times[tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='tc')[0]
    #            trajectory[tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='positions')
                phop     [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='phop')    #  resultsArray[:,tc_index,1]
        #        TAs       [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='TAs')
        #        TBs       [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='TBs')
        #        averageA  [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='averageA')
        #        Prev_StdDev[tc_index]= decoratedFlowRead.read_chunk(frame=tc_index, name='Prev_StdDev')
        #        Next_StdDev[tc_index]= decoratedFlowRead.read_chunk(frame=tc_index, name='Next_StdDev')
        #        average_distAB[tc_index]= decoratedFlowRead.read_chunk(frame=tc_index, name='average_distAB')
        #        softness  [tc_index] = decoratedFlowRead.read_chunk(frame=tc_index, name='softness')
            

                positions = decoratedFlowRead.read_chunk(frame=tc_index, name='positions') # trajectory[tc_index]
                
                descriptor_cutoff = np.max(descriptor_family_parameters[:,0]) - 1e-6
                ## choice of the cutoff: it will be assumed that the max of the 0th column of parmaeters gives a meaningful radius.    

                #positions = Full_positions[selectedAtoms]  # XXX ?
                ############################################################################
                ################ CORE FUNCTION: ############################################
                resultArray, dist_num = module_dist_compute.dist_computer(positions, Lx, expectedDensity, descriptor_cutoff)
                ## we now have all the pair distances from that snapshot. 

                all_atoms_descriptors_array = module_descriptors.descriptor_func_allPositions_oneSnapshot(resultArray, dist_num, descriptor_family_name, descriptor_family_parameters, atomTypes)
                ############################################################################ 
             
                XXX
                TODO
                here the computation I had of the descriptors is most likely qwrong somewhere above... which then gives fucked up, very negative  softness(es)
                
                ## XXX 
                # at this point one may want to save the DESCRIPTORS of the atoms of this frame... that's a lot of data.
                ## we can aslo jjust comopute softness for it. 
                ## XXX
                    
                all_atoms_Softness = np.zeros(Natoms)
                for atom in range(Natoms):
                    X_all_descriptors_array = all_atoms_descriptors_array[atom]
                    
                    X_all_descriptors_array -= descriptor_averages
                    X_all_descriptors_array /= descriptor_stdDevs
                    
                    all_atoms_Softness[atom] = np.dot(a_hyperplane, X_all_descriptors_array) + b_intercept

###                ## TODO: save all_atoms_Softness
###                np.savetxt("softness_t="+str(tc_index)+".dat", all_atoms_Softness)
        
                ff.write_chunk(name='tc'  , data=np.array([tc_index], dtype=np.float32) )
                ff.write_chunk(name='softness', data=np.array( all_atoms_Softness, dtype=np.float32))
                 
                ff.write_chunk(name='phop'     , data=np.array(phop[tc_index], dtype=np.float32))
                ff.end_frame();

print("results written to the file :\n", outputName_softness_R)
#                ### DEBUG : check ###
#                y_calc = np.zeros(Nsamples*2)
#                for i in range(Nsamples*2):
#                    y_calc[i] = np.dot(a_hyperplane,X_all_descriptors_array[i])+b_intercept

#                y_clf = clf.decision_function(X_all_descriptors_array)
#                print(y_clf - y_calc)
#                
#                hard_hist, bins = np.histogram(y_calc[:Nsamples], 100)
#                plt.plot(bins[:-1], hard_hist)
#                
#                soft_hist, bins = np.histogram(y_calc[Nsamples:], 100)
#                plt.plot(bins[:-1], soft_hist)
#                
#                total_hist, bins = np.histogram(y_calc[:], 100)
#                plt.plot(bins[:-1], total_hist)
            

#    decoratedFlowRead.close()
#    print("... file loaded. Use larger every_recAlways if loading is very slow (you may have stupidly set it to 1).")






############################################################################################

########
#using the "table"  pair potential in Hoomd, I can define V and F not from a table, but as functions. 
#And then parameter 'width' says how many discretization points are needed.
#########
#I need to put in the potential term, all the M descriptor functions, and then not just sum them, but compute the distance between
#- this point in d=M dimensions 
#- the plane that was computed once and for all, and is stored in a file or 2, (as a vector + a point in space, to set the reference).
#Of course the parameters used to build the families fo fdescriptors shoudl also be present in a dictionary-like file.

#################################
##### example from the docs:
#def lj(r, rmin, rmax, epsilon, sigma):
#    V = 4 * epsilon * ( (sigma / r)**12 - (sigma / r)**6);
#    F = 4 * epsilon / r * ( 12 * (sigma / r)**12 - 6 * (sigma / r)**6);
#    return (V, F)
#table = pair.table(width=1000)
#table.pair_coeff.set('A', 'A', func=lj, rmin=0.8, rmax=3.0, coeff=dict(epsilon=1.5, sigma=1.0))
#table.pair_coeff.set('A', 'B', func=lj, rmin=0.8, rmax=3.0, coeff=dict(epsilon=2.0, sigma=1.2))
#table.pair_coeff.set('B', 'B', func=lj, rmin=0.8, rmax=3.0, coeff=dict(epsilon=0.5, sigma=1.0))

