import hoomd
from hoomd import md
import sys
import numpy as np
import module_filenamesHandler 
import module_MD_potentials
import gsd.pygsd
import gsd.hoomd

import time

######
## plots: 
#matplotlib.use('Agg') 
import module_importPlotParams
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True          
plt.rcParams.update(module_importPlotParams.fparams())
plt.ion()   ## optional
######


#dt=0.0025
## on GPU, ~1e3 steps are done in ~1sec. 
print 'usage:  \nrun restart.py  (filename.gsd [overriden by "restart.gsd"])   tstepsThermostat   analyzer_adjust_Etot(0 or negative value or number of samples)   tstepsNVE    [glueing_flag:0_no_glueNVT_Fr=XX]    [recPeriod=4e3]    [analyzer_period=4e3]   [tauT=1]  [dt=0.0025]   [thermostat=NVT - or berendsen,brownian,langevin,MB]  [interaction=LJ,WCA,...]  [limit_hours=13.9*24 hours]'
filename = sys.argv[1] #[filename.gsd is overriden by "restart.gsd"] 
tstepsThermostat=int(float(sys.argv[2]))   # 0 is valid
analyzer_adjust_Etot=float(sys.argv[3])     # negative values mans: skip this step.
tstepsNVE=int(float(sys.argv[4]))   # 0 is valid
recPeriod = 4e3
glueing_flag= 'no'
if len(sys.argv)>5:
    glueing_flag = str(sys.argv[5])
    if glueing_flag == 'glueNVT' or glueing_flag == '0' or glueing_flag == 'no' or glueing_flag[:3] == 'Fr=' :
        if glueing_flag[:3] == 'Fr=' :
            Frame_to_restart_from = int(glueing_flag[3:])
            glueing_flag = 'Frame_restart'
    else:
        print 'invalid keyword for glueing_flag. allowed keywords are:  glueNVT , 0 , no , Fr=123456 '
        print 'the last one Fr=X means you create a restart file using Frame #X in hevytrajectories as starting point.\nExit.\n'
        raise SystemExit
    
if len(sys.argv)>6:
    recPeriod = int(float(sys.argv[6]))
analyzer_period = 4e3
if len(sys.argv)>7:
    analyzer_period = int(float(sys.argv[7]))
tauT=1 # float(sys.argv[2])
if len(sys.argv)>8:
    tauT = float(sys.argv[8])
rootname = module_filenamesHandler.get_rootname(filename)
dt=0.0025
if len(sys.argv)>9:
    dt = float(sys.argv[9])
    if dt!=0.0025: 
        if 'dt=' not in filename:
            rootname+="_dt="+str(dt)
        else: 
            dt=float(module_filenamesHandler.filename_parser(filename[:-4], 'dt'))
            #rootname+="_dt="+str(dt)
            #print "accepting, for once, the incoming"
module_filenamesHandler.log_action(rootname, '\n\n'+str(sys.argv)+'\n')
TemperatureGoal = float(module_filenamesHandler.filename_parser(filename[:-4], 'T'))
#TemperatureGoal = 0.23 #45
notice_level=1
thermostat='NVT'
if len(sys.argv)>10:
    thermostat = str(sys.argv[10])
    if thermostat == 'MB':
        notice_level=0
interactionType="LJ"
if len(sys.argv)>11:
    interactionType=sys.argv[11]
    if interactionType == "LJ":
        if 'KA' not in filename:
            print("please label your filenames consistently")
            raise SystemExit
    elif interactionType =="WCA":
        if 'WCA' not in filename:
            print("please label your filenames consistently")
            raise SystemExit
limit_hours=int(24*13.9) # 13.9 days
if len(sys.argv)>12:
    limit_hours=int(float(sys.argv[12]))
#steps_per_tau_Hoomd =  1./XXXXdt
#steps_per_tau_KA    = (1./dt)/ 48**0.5
##steps_per_tau_alpha =((1./dt)/ 48**0.5) * 1e4   # if e.g. tau_alpha = 1e4 in KA units

################################################################################
######## reading the state from a restart file #################################
hoomd.context.initialize(' --notice-level='+str(notice_level)) # --mode=gpu


## glueing_flag == 'Frame_restart' : this is very special, we just read a frame and create a restart file from it, no run performed ##
if glueing_flag == 'Frame_restart':
    if 'HeavyTraj' not in filename:
        print 'error: you are using the mode  glueing_flag=Frame_restart  but you are not providing a "HeavyTrajectory" kind of input file !'
        raise SystemExit
    with open(filename, 'rb') as flow:
        HoomdFlow = gsd.pygsd.GSDFile(flow)
        hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
        Nframes = len(hoomdTraj)
        hoomdTraj.close()
    print 'there are Nframes=', Nframes, ' in the input file'
    if Frame_to_restart_from > Nframes :
        print 'this frame is after the last frame in the input file. Try less.'
        raise SystemExit
    system = hoomd.init.read_gsd(filename=filename, frame = Frame_to_restart_from)
    additional_tag= "-Fr-"+str(Frame_to_restart_from)
    module_filenamesHandler.log_action(rootname, 'Read frame='+str(Frame_to_restart_from)+'  from the HeavyTrajectory and made it into a restart file.\n' )   
    ########## MAJOR CHANGE : ##########
    previous_rootname=rootname
    rootname += additional_tag
    ####################################
    gsd_restart = hoomd.dump.gsd(filename=rootname+"_type=restart.gsd" , group=hoomd.group.all(), period=None, phase=0, truncate=True)  
    gsd_restart.write_restart()
    module_filenamesHandler.log_action(rootname, 'Initial rootname was '+previous_rootname+'. This fork was created by reading the frame='+str(Frame_to_restart_from)+'  from the HeavyTrajectory and making it into a restart file.\n')
    print 'We created the file series: '+rootname+"_type=restart.gsd"
    raise SystemExit    
    ## here we could also continue, but doing so seems more complicated 
    ## than just running again this script with the new restart file as input argument
    ##################### END OF PROGRAM #######################################
#else: 
#    pass

## normal way to restart ##
system = hoomd.init.read_gsd(filename=filename, restart=rootname+'_type=restart.gsd')
restart_period = min(recPeriod, analyzer_period)
gsd_restart = hoomd.dump.gsd(filename=rootname+"_type=restart.gsd" , group=hoomd.group.all(), period=restart_period, phase=0, truncate=True)  

######################################
############# analyzers ##############
######################################

analyzerManyVariables_header     = "# we used (at least initially) dt="+str(dt)+"   See _notes.txt for more comments. "
analyzerManyVariables_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'translational_kinetic_energy', 'rotational_kinetic_energy'] #, 'volume', 'num_particles']
analyzer_NVT_filename                = rootname+"_type=analyzer_"+thermostat+"_tauT="+str(tauT)+".dat"


################################################################################
########## Set up the interactions #############################################
NeighborsListLJ = md.nlist.cell()
if interactionType == "LJ":
    myLjPair = module_MD_potentials.set_interaction_potential_LJ(NeighborsListLJ)    
if interactionType == "WCA":
    myLjPair = module_MD_potentials.set_interaction_potential_WCA(NeighborsListLJ)    
################################################################################
############ Physical integrations #############################################

md.integrate.mode_standard(dt=dt) ##allows Nose-Hoover in NVT, for instance
limit_multiple=int(max(recPeriod, restart_period))
Natoms = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'N')))
hoomd.md.update.zero_momentum(int(1e9/Natoms), phase=0)
#hoomd.md.update.zero_momentum(int(1e6/Natoms), phase=0)
#hoomd.md.update.zero_momentum(10, phase=0)

### other possibilities: 
##md.integrate.berendsen 	Applies the Berendsen thermostat.
##md.integrate.brownian 	Brownian dynamics.
##md.integrate.langevin 	Langevin dynamics.
##md.integrate.mode_standard 	Enables a variety of standard integration methods.


######################################
############# integrators ############
######################################

module_filenamesHandler.log_action(rootname, "Runs (all) starting: ... (currently we are the step  #"+str(hoomd.get_step())+" of the simulation)\n")  ## we put a space in case a run is aborted, to have notes.txt look cleaner.
t0INIT=time.time()

######## NVT ########
if tstepsThermostat > 0 :

    if glueing_flag == 'glueNVT' :
        curStep = hoomd.get_step()
        gsd_HeavyTrajectory = hoomd.dump.gsd(filename=rootname+"_type=HeavyTraj"+"_cst="+str(int(curStep/1e9))+"e9_rP="+str(int(recPeriod))+".gsd", group=hoomd.group.all(), period=recPeriod, phase=0, overwrite=False, truncate=False, static=['attribute', 'topology'])
        ## we take out "momentum" from the static qties, so it will be recorded at each frame
    ## end of flag ## 

    ## hack to anneal a simu ##
    #TemperatureGoal=5.0
    #TemperatureGoal=1.0
    ## end of hack ##
    
    if thermostat == 'berendsen' :
        integrator_nvt = md.integrate.berendsen(group=hoomd.group.all(), kT=TemperatureGoal, tau=tauT)
        module_filenamesHandler.log_action(rootname, "Starting a Berendsen (velocity rescaling) run for a number  "+str(tstepsThermostat)+" == 10**"+str(np.round(np.log10(tstepsThermostat),1))+"  of time steps, with dt="+str(dt)+" and tauT="+str(tauT)+" : ... (run starts)")
    elif thermostat == 'brownian' :
        integrator_nvt = md.integrate.brownian(group=hoomd.group.all(), kT=TemperatureGoal, seed=42, dscale=tauT)
        module_filenamesHandler.log_action(rootname, "Starting a Brownian (overdamped) run for a number  "+str(tstepsThermostat)+" == 10**"+str(np.round(np.log10(tstepsThermostat),1))+"  of time steps, with dt="+str(dt)+" and tauT="+str(tauT)+" : ... (run starts)")
    elif thermostat == 'MB' :
        stepsTauT = int(tauT/dt)
        integrator_nvt = md.integrate.nve(group=hoomd.group.all())

        module_filenamesHandler.log_action(rootname, "Starting a MaxBolt (MB+velocity rescaling) run for a number  "+str(tstepsThermostat)+" == 10**"+str(np.round(np.log10(tstepsThermostat),1))+"  of time steps, with dt="+str(dt)+" and tauT="+str(tauT)+" : ... (run starts)")
    elif thermostat == 'NVT' :
        integrator_nvt = md.integrate.nvt(group=hoomd.group.all(), kT=TemperatureGoal, tau=tauT)
        module_filenamesHandler.log_action(rootname, "Starting a NVT run for a number   "+str(tstepsThermostat)+" == 10**"+str(np.round(np.log10(tstepsThermostat),1))+"  of time steps, with dt="+str(dt)+" and tauT="+str(tauT)+" : ... (run starts)")

    analyzerManyVariables_NVT = hoomd.analyze.log(filename=analyzer_NVT_filename, \
            quantities=analyzerManyVariables_quantities, period=analyzer_period, \
            header_prefix = analyzerManyVariables_header+' this is the NVT part, as said in the title' , overwrite=False, phase=0)

    if thermostat == 'MB' :
        t0=time.time()
        for iterations in range(0,int(tstepsThermostat/stepsTauT)):
            
            ##reset velocities
            snap = system.take_snapshot()
            # each component is a gaussian of variance sigma, that's it.
            vel = np.random.normal(0,TemperatureGoal**0.5, (Natoms,3))  
#            print '\nVelcities were rescaled. by ',TemperatureGoal**0.5/np.std(vel), '\n'
            vel *= TemperatureGoal**0.5/np.std(vel)
            snap.particles.velocity[:] = vel
            system.restore_snapshot(snap)            
            hoomd.run(stepsTauT, quiet=False)
        t1=time.time()
    else: 
        #hoomd.run_upto(tstepsThermostat, quiet=False) 
        t0=time.time()
        hoomd.run(tstepsThermostat, quiet=False) 
        t1=time.time()
    module_filenamesHandler.log_action(rootname, " (run ends)... "+thermostat+" thermostat  completed.\n(it took  "+str(int(t1-t0))+" seconds, i.e.  "+str(np.round((t1-t0)/3600,1))+"  hours)\n")

    analyzerManyVariables_NVT.disable()
    integrator_nvt.disable()


if analyzer_adjust_Etot != 0 : 
#####  Daniele Coslovitch trick: when unolugging the thermostat, adjust the Ec such that E_tot is equal to the target value:
### <Etot_before> = Etot_target. Etot = Ec + Ep. Recale Ec to (Etot_target-Ep), i.e. rescale velocities by:
### Vnew = v_now * (<Etot_before>-Ep_now)/Ec_now. 
###- if you do not adjust Etot, the system will adjust its Ep and Ec to have equipartition such that the actual T of your new simulation will be slightly off the target Ec (or T). You can als accept that and simply mention in your paper "my temperature was 0.43483" instead of "0.43". For low T it seems nicer to adjust things.

    ### computing the target value of Etot:  ###
    if analyzer_adjust_Etot > 0 :
# TODO: do not use more than the very prevous NVT run, except if it was steps=0 (ovrride): otherwise alternating NVT and NVE may causes gaps in Etot between the 2 NVT runs that are not consecutive .. and then the rescaling is bad.
        tout = np.loadtxt(analyzer_NVT_filename)
        if len(tout) > analyzer_adjust_Etot :
            tout = tout[-int(analyzer_adjust_Etot):, :]
        else: 
            module_filenamesHandler.log_action(rootname,"\nError: less data available than asked. (you have "+str(len(tout))+" and you requested "+str(int(analyzer_adjust_Etot))+" \n")
            raise SystemExit
        step,T,p,Ep,Ec,M = tout[:,0], tout[:,1], tout[:,2], tout[:,3], tout[:,4], tout[:,5]
        stdDevTemperature = np.std(T)
        if stdDevTemperature > 0.1 :
            module_filenamesHandler.log_action(rootname,"\n\nWARNING: your temperature fluctates a lot in the NVT regime you are using to adjust E_tot.  StdDev)T)="+str( stdDevTemperature)+'\n\n')
#            del tout
        Etot_target = np.mean(Ep+Ec)
        module_filenamesHandler.log_action(rootname,"\nComputed <Etot> = "+str( Etot_target)+ " from a number N="+str( analyzer_adjust_Etot)+ " data points.")
    else: 
        Etot_target = analyzer_adjust_Etot
    ## now we have Etot_target

    ## measuring Ec, Ep to be able to adjust Ec ##
    integrator_adjust = md.integrate.nve(group=hoomd.group.all())
    analyzerManyVariables_NVT = hoomd.analyze.log(filename='adjust.dat', \
        quantities=analyzerManyVariables_quantities, period=1, \
        header_prefix = analyzerManyVariables_header+' this is the NVT part, as said in the title' , overwrite=False, phase=0)
    hoomd.run(1, quiet=True)
    TempNOW, EpNOW, EcNOW = analyzerManyVariables_NVT.query('temperature') , analyzerManyVariables_NVT.query('potential_energy'),  analyzerManyVariables_NVT.query('kinetic_energy')
    
    longWhileLoop=0
    while abs((EcNOW+EpNOW)/Etot_target-1)> 1e-4 :  ## relaive error should be small.
        longWhileLoop+=1
        if longWhileLoop > 100:
            module_filenamesHandler.log_action(rootname,'\nAdjust is not convergng after 100 idteration... lets leave ! \n') 
            break
#        print("\n\n")
#        print ("(EcNOW+EpNOW)/Etot_target-1 = ", (EcNOW+EpNOW)/Etot_target-1 )
#        print("TempNOW, EpNOW, EcNOW", TempNOW, EpNOW, EcNOW)

        snap = system.take_snapshot()
        vel = snap.particles.velocity[:]
#        EcManual = 0.5*1*np.sum(vel**2)
#        print("EcManual=", EcManual,  "  but EcNOW=",EcNOW)
######        EcNOW=EcManual

    #        Natoms = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'N')))
    #        print(rootname,EcNOW , " ==??== ",  3./2. * (Natoms-1) * TempNOW, " ==??== ", 3 * 0.5 * np.var(vel) *(Natoms)  ) ## 1/2 m v^2 becomes 3/2 because each component is summed independently in 3D . ##  understand the origin of the -1 in Natoms-1
        vel *=  ( (Etot_target - EpNOW) / (EcNOW) )**0.5
        if np.sum( (vel>-1000000000000) ) == np.size(vel) : ## check there is no nan, inf and so on
            module_filenamesHandler.log_action(rootname,'\nwe rescale velocities by a factor = '+str( ( (Etot_target - EpNOW) / (EcNOW) )**0.5 )+ ' ') 
            snap.particles.velocity[:] = vel
            system.restore_snapshot(snap)
        else:
            module_filenamesHandler.log_action(rootname,'\nwe FAIL to rescale velocities by a factor = '+str( ( (Etot_target - EpNOW) / (EcNOW) )**0.5 )+ ',  because we have some nan..') 
            pass
#        hoomd.run(limit_multiple, quiet=True)
        hoomd.run(1, quiet=True)
        TempNOW, EpNOW, EcNOW = analyzerManyVariables_NVT.query('temperature') , analyzerManyVariables_NVT.query('potential_energy'),  analyzerManyVariables_NVT.query('kinetic_energy')
#        print("TempNOW, EpNOW, EcNOW", TempNOW, EpNOW, EcNOW)
#        snap = system.take_snapshot()
#        vel = snap.particles.velocity[:]
#        EcManual = 0.5*1*np.sum(vel**2)
#        print("EcManual=", EcManual,  "  but EcNOW=",EcNOW)

    hoomd.run(recPeriod-1-longWhileLoop)
    analyzerManyVariables_NVT.disable()
    module_filenamesHandler.log_action(rootname,"\nNow we have :  (EcNOW+EpNOW)/Etot_target-1 = " +str( (EcNOW+EpNOW)/Etot_target-1 ) +'\n\n')
    integrator_adjust.disable()


######## NVE ########
if tstepsNVE > 0 :
    curStep = hoomd.get_step()
    gsd_trajectory = hoomd.dump.gsd(filename=rootname+"_type=traj"+"_cst="+str(int(curStep/1e9))+"e9_rP="+str(int(recPeriod))+".gsd", group=hoomd.group.all(), period=recPeriod, phase=0, overwrite=False)

    analyzer_NVE_filename = rootname+"_type=analyzer_NVE.dat"
    analyzerManyVariables_NVE = hoomd.analyze.log(filename=analyzer_NVE_filename, \
            quantities=analyzerManyVariables_quantities, period=analyzer_period, \
            header_prefix = analyzerManyVariables_header+' this is the NVE part, as said in the title' , overwrite=False, phase=0)

    integrator_nve = md.integrate.nve(group=hoomd.group.all())
    module_filenamesHandler.log_action(rootname, "Starting a NVE run for a number   "+str(tstepsNVE)+" == 10**"+str(np.round(np.log10(tstepsNVE),1))+"  of time steps, with dt="+str(dt)+"(at constand energy) : ... (run starts)")
    #hoomd.run_upto(tstepsThermostat+tstepsNVE, profile=False, limit_hours=limit_hours, limit_multiple=max(recPeriod, restart_period), quiet=False) 
    t0=time.time()
    hoomd.run(tstepsNVE, quiet=False, profile=False, limit_hours=limit_hours, limit_multiple=limit_multiple)
    t1=time.time()
    module_filenamesHandler.log_action(rootname, " (run ends)... NVE run completed.\n(it took  "+str(int(t1-t0))+" seconds, i.e.  "+str(np.round((t1-t0)/3600,1))+"  hours)\n")
    analyzerManyVariables_NVE.disable()
    integrator_nve.disable()

t1FIN=time.time()
module_filenamesHandler.log_action(rootname, "\n(all the run(s) end)... all runs finished. (currently we are the step  #"+str(hoomd.get_step())+" of the simulation)\n(Together, they took  "+str(int(t1FIN-t0INIT))+" seconds, i.e.  "+str(np.round((t1FIN-t0INIT)/3600,1))+"  hours)\n\n")  ## we put a space in case a run is aborted, to have notes.txt look cleaner.





if tstepsNVE > 0 :
#try : 
    toplot = np.loadtxt(analyzer_NVE_filename)
    plt.figure(1,[20,6])
    steps = toplot[:,0]*dt
    T = toplot[:,1]
    plt.plot(steps, T, ls= ' ', marker='+', label='$T$, NVE')
    plt.legend(loc='upper right')
    plt.xlabel('$t$ (LJ units)')
    plt.ylabel(r'$T$')
    outName=rootname + "_ANALYZ=T-NVE_analyzer"
    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
#    plt.close()

    plt.figure(2,[20,6])
    toplot = np.loadtxt(analyzer_NVE_filename)
    steps = toplot[:,0]*dt
    Etot = toplot[:,3]+toplot[:,4]
    plt.plot(steps, Etot, ls= ' ', marker='x', label='$E_{tot}$, NVE')
    plt.legend(loc='upper right')
    plt.xlabel('$t$ (LJ units)')
    plt.ylabel(r'$E_{tot}$')
    outName=rootname + "_ANALYZ=Etot-NVE_analyzer"
    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
#    plt.close()
#except:
#    print("fail to plot NVE-related E and T")
#    pass

try: 
#if tstepsThermostat > 0 :
    plt.figure(1,[20,6])
    toplot = np.loadtxt(analyzer_NVT_filename)
    steps = toplot[:,0]*dt
    T = toplot[:,1]
    plt.plot(steps, T, ls= ' ', marker='+', label='$T$, '+thermostat+r' at $\tau_T='+str(tauT)+'$')
    plt.legend(loc='upper right')
    plt.xlabel('$t$ (LJ units)')
    plt.ylabel(r'$T$')
    outName=rootname + "_ANALYZ=T-NVT_analyzer"
    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
    plt.close()

    plt.figure(2,[20,6])
    toplot = np.loadtxt(analyzer_NVT_filename)
    steps = toplot[:,0]*dt
    Etot = toplot[:,3]+toplot[:,4]
    plt.plot(steps, Etot, ls= ' ', marker='x', label='$E_{tot}$, '+thermostat+r' at $\tau_T='+str(tauT)+'$')
    plt.legend(loc='upper right')
    plt.xlabel('$t$ (LJ units)')
    plt.ylabel(r'$E_{tot}$')
    outName=rootname + "_ANALYZ=Etot-NVT_analyzer"
    module_importPlotParams.savefig_perso(outName+module_importPlotParams.fsaveFigFormat())
    plt.close()
except:
    print("fail to plot thermostat E and T")
    pass
    

#################################################################################
#################################################################################
#################################### The End ####################################

#gsd_restart.write_restart()


#            rescaled_cheat = md.integrate.berendsen(group=hoomd.group.all(), kT=TemperatureGoal, tau=dt*1.0)
#            rescaled_cheat.enable()
#            hoomd.run(100, quiet=False)
#            rescaled_cheat.disable()
#            integrator_nvt = md.integrate.nve(group=hoomd.group.all())
#            integrator_nvt.enable()

