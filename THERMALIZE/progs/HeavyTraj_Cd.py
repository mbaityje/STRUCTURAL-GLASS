import hoomd
from hoomd import md
from lib import module_potentials as pot
import numpy as np
# from lib import module_measurements as med


context=hoomd.context.initialize()


#PARAMETERS
T=2.0
Natoms=1080
maindir='../../OUTPUT/T'+str(T)+'/N'+str(Natoms)+'/'
samples=range(10)
nsamples=len(samples)
EPS=1e-8


#READ THERMALIZED SYSTEM (for copying particle types and creating snapshots) AND POTENTIAL
sam=0
system = hoomd.init.read_gsd(filename=maindir+'/S'+str(sam)+'/thermalized.gsd')
L=system.box.Lx
pair=pot.LJ(md.nlist.cell(), type='KAshort')


#READ HEAVY TRAJECTORY

#Make sure that L and dt are consistent
for isam in range(nsamples):
	sam=samples[isam]
	name=maindir+'/S'+str(sam)+'/heavyTraj/L.txt'
	f=open(name,'rt')
	Ltemp=np.float64(f.readline())
	if 0==isam:
		dt=np.float64(f.readline())
	else:
		dttemp=np.float64(f.readline())
		if np.abs(dt-dttemp)>EPS: raise ValueError('isam %d has an incompatible dt (%.9g) with sample 0 (%.9g)'%(isam,dttemp,dt))
	if np.abs(L - Ltemp)>EPS: raise ValueError('isam %d has an incompatible L  (%.9g) with sample 0 (%.9g)'%(isam, Ltemp, L))
	f.close()



#Read positions, velocities and accelerations of all samples
all_times=[]
all_pos=[]
all_vel=[]
all_acc=[]
ntimes=np.zeros(nsamples,dtype=int)
for isam in range(nsamples):
	sampledir=maindir+'/S'+str(isam)+'/heavyTraj/'
	#Read list of times
	timesName=sampledir+'times.txt'
	timelist=np.loadtxt(timesName)
	ntimes[isam]=len(timelist)
	print('There are ',ntimes[isam],' configurations in total')

	
	#Read accelerations
	pos = []
	vel = []
	acc = []
	fpos=open(sampledir+'/pos.npy','rb')
	fvel=open(sampledir+'/vel.npy','rb')
	facc=open(sampledir+'/acc.npy','rb')
		
	for i in range(ntimes[isam]):
		pos.append(np.load(fpos))
		vel.append(np.load(fvel))
		acc.append(np.load(facc))
	fpos.close()
	fvel.close()
	facc.close()

	
	#Relevant quantities
	if 0==isam:
		Natoms=len(pos[0])
	elif Natoms!=len(pos[0]): 
		raise ValueError('isam %d has an inconsistent Natoms=%d, inconsistent with sample 0 that has %d'%(isam,len(pos[0]),Natoms))
	initialPositions=pos[0]
	all_times.append(timelist)
	all_pos.append(pos)
	all_vel.append(vel)
	all_acc.append(acc)

del pos,vel,acc,timelist
	
all_pos=np.array(all_pos)
all_vel=np.array(all_vel)
all_acc=np.array(all_acc)





#Waiting times
#snumber of trajectory starting times for each sample
ntw=np.zeros(nsamples,dtype=np.int)
twlist=[]
for isam in range(nsamples):
	twlist.append(np.unique(all_times[isam][:,0]))
	ntw[isam]=len(twlist[isam])   
print("ntw: ",ntw)



# Sort Observables
nt=np.zeros(nsamples,dtype=np.int)
tlist  =[]
poslist=[]
vellist=[]
acclist=[]
print(np.shape(all_pos))
for isam in range(nsamples):    
	tlist_sam  =[]
	poslist_sam=[]
	vellist_sam=[]
	acclist_sam=[]
	nt[isam]=0
	for itw in range(ntw[isam]):
		tw=twlist[isam][itw]
		this_tw=np.where(all_times[isam][:,0]==tw)
		first=this_tw[0][0]
		last=this_tw[0][-1]+1
		if last-first>nt[isam]:
			nt[isam]=last-first
			tlist_sam=np.array(all_times[isam][first:last][:,1]-all_times[isam][first][1], dtype=np.int64)
		poslist_sam.append(all_pos[isam][first:last])
		vellist_sam.append(all_vel[isam][first:last])
		acclist_sam.append(all_acc[isam][first:last])
	
	#Remove incomplete entries
	for itw in np.arange(ntw[isam]-1,-1,-1):
		if len(poslist_sam[itw]) != nt[isam]:
			assert(len(vellist_sam[itw]) == len(poslist_sam[itw]))
			assert(len(acclist_sam[itw]) == len(poslist_sam[itw]))
			del poslist_sam[itw]
			del vellist_sam[itw]
			del acclist_sam[itw]
			twlist_sam=np.delete(twlist_sam,itw)

	#poslist, vellist and acclist have the shape: [itw, it, particle, component]

	print('ntw:',ntw[isam],end='  ')
	ntw[isam]=len(twlist[isam])
	print('ntw:',ntw[isam])

	tlist.append(tlist_sam)
	poslist.append(poslist_sam)
	vellist.append(vellist_sam)
	acclist.append(acclist_sam)
del poslist_sam,vellist_sam,acclist_sam,tlist_sam

#The following objects should be of size [nsamples, ntw, nt, Natoms, Dim]
print('ntw must be the same for all samples: ',ntw)
print('nt  must be the same for all samples: ',nt)
print(np.shape(poslist))
print(np.shape(vellist))
print(np.shape(acclist))









###################
# CALCULATE Cd(t) #
###################
snapA=system.take_snapshot()
snapB=system.take_snapshot()

n1=1; n2=3; n3=4

obs=[{} for isam in range(nsamples+1)]
Cd = np.zeros((nsamples, max(ntw), max(nt)), dtype=np.float64)
for isam in range(nsamples):
	for itw in range(ntw[isam]):
		initialPositions[:]=poslist[isam][itw][0][:]
		snapA.particles.position[:] = initialPositions[:]
		for iframe in range(0, nt[isam]):
			snapB.particles.position[:] = poslist[isam][itw][iframe][:]
			Cd[isam][itw][iframe]=pair.Cd(snapA=snapA,snapB=snapB,beta=1./T)
			print('isam: %d, itw: %d, it: %d, Cd= %g\n'%(isam,itw,iframe,Cd[isam][itw][iframe]))

	obs[isam]={
        'Cd':{'mean': np.mean(Cd[isam],axis=0), 'err': sem(Cd[isam], axis=0)},
    }





