#!/usr/bin/env python
# python IdentifyMetabasins.py ~/STRUCTURAL-GLASS/OUTPUT/T0.6/N65/shift/S0/chunksIS/elistIS.txt

import sys
import numpy as np, pandas as pd
import argparse, itertools
import lib.module_measurements as med
from matplotlib import pyplot as plt

# READ COMMAND-LINE ARGUMENTS
parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
parser.add_argument('dirname', help='Name of the file with the IS energy')
parser.add_argument('--dt', type=float ,default=0.0025, help='Time step (only for figures)')
parser.add_argument('--thres', type=float ,default=1e-4, help='Threshold for metabasin energy precision')
parser.add_argument('--maxbasins', type=int ,default=10, help='Measure the deepest maxbasins MetaBasins')
parser.add_argument('--mintau', type=int ,default=0, help='Smallest basin persistence time. Default is zero (i.e. basin immediately abandoned)')
parser.add_argument('-L','--L', type=np.float64, default=3.78364777565, help='Size of the box. Default is 3.78364777565, the value for the N=65, rho=1.2 systems')
parser.add_argument('--showplots', action='store_true', help='Flag to show plots')
parser.add_argument('--noLJunits', action='store_true', help='Flag to show time step in number of MD iterations instead of LJ time')
args = parser.parse_args()


# READ LIST OF INHERENT STRUCTURES
eISname=args.dirname+'/elistIS.txt'
df=pd.read_csv(eISname, sep=' ',header=None, names=['time','Eis'], skiprows=1, dtype={'time':np.float64,'Eis':np.float64})
nlines=len(df['time'])
tmax=df['time'].iloc[nlines-1]


# READ LIST OF RIDGES
eRidgeName=args.dirname+'/elistRidge.txt'
dfr=pd.read_csv(eRidgeName, sep=' ',header=None, names=['time','Eridge','Eante','Epost'], skiprows=1, dtype={'time':np.float64,'Eridge':np.float64,'Eante':np.float64,'Epost':np.float64})
# nlinesr=len(dfr['time'])
# tmaxr=dfr['time'].iloc[nlines-1]


def t_from_i(i):
	return np.int64(df['time'].iloc[i])

def appendToIntervals(intervals, ia, ib, ini, fin):
	ea=dftemp['Eis'].iloc[ia]
	eb=dftemp['Eis'].iloc[ib]
	ebasin=dftemp['Eis'].iloc[ini]

	# In the following, we treat the case in which the initial or final configuration
	# of the interval has a different (higher) energy. In that case, we reduce the size
	# of the interval so that beginning and end of the interval have the same energy

	if ( np.abs( ea - ebasin) > THRES ): # Piscio l'intervallo se l'inizio non ha l'energia giusta
		for i in range(ia,ib):
			if ( np.abs( dftemp['Eis'].iloc[i] - ebasin) < THRES ):
				ia=i
				break

	if ( np.abs( eb - ebasin) > THRES ): # Piscio l'intervallo se l'inizio non ha l'energia giusta
		for i in range(ib,ia-1,-1):
			if ( np.abs( dftemp['Eis'].iloc[i] - ebasin) < THRES ):
				ib=i
				break
		return

	intervals.append( (ia, ib ) )
	return

dftemp=df.copy(deep=True)

basins=[]

MAXBASINS=args.maxbasins
THRES=args.thres
ANALYZED=1
for i in range(MAXBASINS):
	minimo=dftemp['Eis'].min()
	if minimo>0: 
		break
	imins=dftemp.index[ np.abs(dftemp['Eis']-minimo)<THRES].tolist()

	ini=int(imins[ 0])
	fin=int(imins[-1])
	t_ini=np.int64(df['time'].iloc[ini]) if fin<nlines else np.int64(df['time'].iloc[nlines-1])
	t_fin=np.int64(df['time'].iloc[fin]) if fin<nlines else np.int64(df['time'].iloc[nlines-1])+1
	
	if np.abs(minimo-dftemp['Eis'].iloc[ini])>THRES:
		raise ValueError('minimo: %g   elem: %g   |diff|: %g'%(minimo,dftemp['Eis'].iloc[ini], np.abs(minimo-dftemp['Eis'].iloc[ini])))
	if np.abs(minimo-dftemp['Eis'].iloc[fin])>THRES:
		raise ValueError('%g  %g'%(minimo,dftemp['Eis'].iloc[fin]))

	if t_fin-t_ini>=args.mintau: #Exclude MBs smaller than args.mintau

		# If there is no deeper basin (i.e. a basin that has already been analyzed) 
		# between those values, we found the basin
		analyzed=np.where(np.abs(dftemp['Eis'].iloc[ini:fin+1] - ANALYZED) <1e-10)[0]
		analyzed=analyzed+ini
		if len(analyzed ) == 0:
			basins.append((minimo, t_ini, t_fin))
			dftemp['Eis'].iloc[ini:fin+1]=ANALYZED

		#If there is a previously-found (deeper) basin between those values, we must consider the subintervals:
		else: 
			intervals=[]
			appendToIntervals(intervals, ini, analyzed[0]-1, ini=ini, fin=fin)
			for i in range(len(analyzed)-1):
				if analyzed[i+1]-analyzed[i]>1:
					appendToIntervals(intervals, analyzed[i], analyzed[i+1], ini=ini, fin=fin )
			appendToIntervals(intervals, analyzed[-1]+1, fin, ini=ini, fin=fin)

			for interv in intervals:
				basins.append((minimo, t_from_i(interv[0]), t_from_i(interv[1])) )
				dftemp['Eis'].iloc[interv[0]:interv[1]+1]=ANALYZED
	else:
		dftemp['Eis'].iloc[ini:fin+1]=ANALYZED
print('Targeted {} MetaBasins'.format(i))

del dftemp
basins=np.array(basins)
basins=basins[basins[:,1].argsort()]


#Proviamo a identificare ste metabasin ridges
MBridgeTimes=basins[:-1][:,2]+0.5
MBridges=[]
for i in range(len(dfr['time'].values[:])):
	if dfr['time'].values[i] in MBridgeTimes:
		MBridges.append([dfr['time'].values[i],dfr['Eridge'].values[i]])
MBridges=np.array(MBridges)




if args.showplots:
	import pylab as pl
	from matplotlib import collections as mc
	LJunits=True

	fig, ax = pl.subplots()

	if args.noLJunits==False:
		lines = [[(row[1]*args.dt,row[0]),(row[2]*args.dt,row[0])] for row in basins]
		lc = mc.LineCollection(lines, linewidths=4, color='black')
		ax.plot(df ['time']*args.dt,df ['Eis'], 'o-',zorder=0)
		ax.plot(dfr['time']*args.dt,dfr['Eridge'], 'o',zorder=0)
		ax.plot(MBridges[:,0]*args.dt,MBridges[:,1], 'x',zorder=0, color='black')
	else:
		lines = [[(row[1],row[0]),(row[2],row[0])] for row in basins]
		lc = mc.LineCollection(lines, linewidths=4, color='black')
		ax.plot(df ['time'], df['Eis'   ], 'o-',zorder=0)
		ax.plot(dfr['time'],dfr['Eridge'], 'o' ,zorder=0)
		ax.plot(MBridges[:,0],MBridges[:,1], 'x',zorder=0, color='black')
	ax.set_xlabel(r'$t$')
	ax.set_ylabel(r'$E_\mathrm{IS}$')
	ax.add_collection(lc)
	ax.autoscale()
	ax.margins(0.1)
	# plt.ylim((-406.5,-395))
	plt.savefig('Eis.png')
	plt.show()



#########################################################################
#                                                                       #
# SECOND PART OF THE SCRIPT: READ CONFIGURATIONS AND CALCULATE OVERLAPS #
#                                                                       #
#########################################################################
import gsd.pygsd
import gsd.hoomd

def FindtChunk(filename):
	''' Finds the number of time steps and frames in the given chunk'''
	with open(filename, 'rb') as flow:
		HoomdFlow = gsd.pygsd.GSDFile(flow)
		hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
		Nframes = len(hoomdTraj)
		s0    = hoomdTraj.read_frame(0)     # Snapshot of initial configuration (frame zero)
		slast = hoomdTraj.read_frame(Nframes-1) # Snapshot of   final configuration

	if s0.particles.N != 65:
		raise ValueError('trajChunk0.gsd has different from Natoms=65. In this project I am doing only N=65')

	return slast.configuration.step, Nframes


#
# Fore every basin, find one instance of every IS, and report the time at which it occurs
#

# Find times of different Inherent Structures (IS)

allbtimes=[] #All the times at which we can find a different IS

skippedbasins=1 #Skip the first and last skippedbasins basins because they are too close to the border of the run
for ib in range(skippedbasins,len(basins)-skippedbasins): #for each MB, find all its different ISs
	bindices = df.loc[ (df['time']>=basins[ib,1]) & (df['time']<=basins[ib,2]) ].index # Indices of the ISs in the metabasin
	beners   = df.iloc[bindices]['Eis'].round(decimals=6).drop_duplicates()            # Each different energy in the metabasin
	btimes   = df.iloc[beners.index]['time'].values.astype(int)                        # The times at which each energy occurs
	allbtimes.extend(btimes)
allbtimes=set(allbtimes)

# Read all IS configurations

chunkTime, lenChunk=FindtChunk(args.dirname+'/trajChunk0.gsd')
nChunks=int(tmax/lenChunk)+1
ib=skippedbasins # Index of basin starts
confbasins={'ibasin':[], 'times':[], 'confs':[]} # index of metabasin, times of its ISs, confs of its ISs
for ichunk in range(nChunks):
	trajname=args.dirname+"/trajIS"+str(ichunk)+".npy"
	in_trajIS=open(trajname, "rb")
	lastTime=chunkTime+lenChunk*ichunk
	for i in range(lenChunk): #This loop should break at i<<maxt, because trajIS is the result of bisection (does not report every single time step as does trajChunk)
		tempo, config=np.load(in_trajIS)

		if tempo == basins[ib][1]:
			assert(tempo in allbtimes)
			confbasins['ibasin'].append(ib-skippedbasins)
			confbasins['times'  ].append([tempo])
			confbasins['confs'  ].append([config])
		elif tempo in allbtimes:
			confbasins['times'  ][-1].append(tempo)
			confbasins['confs'  ][-1].append(config)

		if tempo>=basins[ib][2]:  #If the basin ended, go to the next one
			ib+=1
			if ib==len(basins)-skippedbasins: #If we are beyond the last accounted basin there is no need to keep reading useless data
				break

		if tempo==lastTime: 
			break


#
# Calculate intrabasin and extrabasin overlaps
#
nbasins=len(confbasins['ibasin'])
confbasins['qintra']=np.full(nbasins, None)
confbasins['qextra']=[[None]]*nbasins # Cannot initialize as numpy array, bc every element has different number of subelements
for ib in range(len(confbasins['ibasin'])):

	#Intrabasin overlap
	nconfs_i=len(confbasins['times'][ib])
	qin=0.0
	for ic in range(nconfs_i):
		for jc in range(ic+1,nconfs_i):
			qin += med.OverlapPos(confbasins['confs'][ib][ic], confbasins['confs'][ib][jc], args.L)
	qin /= 0.5*(nconfs_i-1)*nconfs_i if nconfs_i>1 else np.nan
	confbasins['qintra'][ib]=qin

	#Extrabasin overlap with the subsequent basins
	for jb in range(ib+1, nbasins):
		qex=0.0
		nconfs_j=len(confbasins['times'][jb])
		for ic in range(nconfs_i):
			for jc in range(nconfs_j):
				qex += med.OverlapPos(confbasins['confs'][ib][ic], confbasins['confs'][jb][jc], args.L)
		qex/=(nconfs_i*nconfs_j)
		if jb==ib+1:
			confbasins['qextra'][ib] = [qex]
		else:
			confbasins['qextra'][ib].append(qex)




#####################
#                   #
# FINAL OBSERVABLES #
#                   #
#####################

# Metabasin time
tau = np.array([basins[ib,2]-basins[ib,1] for ib in range(skippedbasins, len(basins)-skippedbasins)], dtype=int)

# Intrabasin overlap
qintra = confbasins['qintra'][ np.where(confbasins['qintra']>0)[0] ][0] #The np.where is because single-conf MBs don't have a self-overlap and are set to nan

# Extrabasin overlap with the following basins (last basin is removed because nothing follows it)
# There are nbasins-1 possible distances, and nbasins-1 basins from which to measure a qextra with subsequent MBs
qdist=np.array([np.mean([confbasins['qextra'][ib][dist] for ib in range(nbasins-1-dist)] ) for dist in range(nbasins-1)])

# Extrabasin overlaps. Can't average qevol, otherwise small distances get overrepresented
qextralist = np.array(list(itertools.chain.from_iterable(confbasins['qextra'])))[:-1] #Last element is None, so we exclude it
qextra = qextralist.mean()



#####################
#                   #
#       OUTPUT      #
#                   #
#####################

# Save tau
np.savetxt(args.dirname+'tauMB.txt',tau, fmt="%d")

# Save q
overlaps=list(zip(qdist, np.full(len(qdist),qintra), np.full(len(qdist),qextra)))
np.savetxt(args.dirname+'qMB.txt', overlaps, fmt="%.14g %.14g %.14g")

# Save ridges
np.savetxt(args.dirname+'EridgeMB.txt',MBridges, fmt="%d")

# Save Metabasin Energies
np.savetxt(args.dirname+'EMB.txt',basins[:,0], fmt="%d")

