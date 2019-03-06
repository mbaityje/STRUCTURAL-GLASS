#!/usr/bin/env python
# python IdentifyMetabasins.py ~/STRUCTURAL-GLASS/OUTPUT/T0.6/N65/shift/S0/chunksIS/elistIS.txt

import sys
import numpy as np, pandas as pd
import argparse
import lib.module_measurements as med
from matplotlib import pyplot as plt

# READ COMMAND-LINE ARGUMENTS
parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
parser.add_argument('eISname', help='Name of the file with the IS energy')
parser.add_argument('--dt', default=0.0025, help='Time step (only for figures)')
args = parser.parse_args()


df=pd.read_csv(args.eISname, sep=' ',header=None, names=['time','Eis'], skiprows=1, dtype={'time':np.float64,'Eridge':np.float32,'Eante':np.float32,'Epost':np.float32,})
print(df.tail())
nlines=len(df['time'])


def t_from_i(i):
	return np.int64(df['time'].iloc[i])

dftemp=df.copy(deep=True)

basins=[]

MAXBASINS=100
THRES=1e-4
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
		raise ValueError('minimo: %g   elem: %g   |diff|: %g'%(minimo,dftemp['Eis'].iloc[ini],np.abs(minimo-dftemp['Eis'].iloc[ini])))
	if np.abs(minimo-dftemp['Eis'].iloc[fin])>THRES:
		raise ValueError('%g  %g'%(minimo,dftemp['Eis'].iloc[fin]))
	if t_fin-t_ini>0:
		# If there is no deeper basin (i.e. a basin that has already been analyzed) 
		# between those values, we found the basin
		analyzed=np.where(np.abs(dftemp['Eis'].iloc[ini:fin+1] - ANALYZED) <1e-10)[0]
		analyzed=analyzed+ini
		if len(analyzed ) == 0:
			basins.append((minimo, t_ini, t_fin))
			#If there is previously found (deeper) between those values, we must consider the subintervals
		else:
			intervals=[(t_ini, t_from_i(analyzed[0]))]
			for i in range(len(analyzed)-1):
				if analyzed[i+1]-analyzed[i]>1:
					intervals.append((t_from_i(analyzed[i]), t_from_i(analyzed[i+1]) ))
			intervals.append((t_from_i(analyzed[-1]), t_fin))

			for interv in intervals:
				basins.append((minimo, interv[0], interv[1]))
		 
	dftemp['Eis'].iloc[ini:fin+1]=ANALYZED

# del dftemp
basins=np.array(basins)
basins=basins[basins[:,1].argsort()]


import pylab as pl
from matplotlib import collections as mc

lines = [[(row[1]*args.dt,row[0]),(row[2]*args.dt,row[0])] for row in basins]
fig, ax = pl.subplots()
lc = mc.LineCollection(lines, linewidths=4, color='black')
ax.plot(df['time']*args.dt,df['Eis'],zorder=0)
ax.add_collection(lc)
ax.autoscale()
ax.margins(0.1)
# plt.ylim((-406.5,-395))
plt.show()
