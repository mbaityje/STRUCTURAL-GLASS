#!/usr/bin/env python
#
# Module for linear and logarithmically binned histograms
# 
# 

import numpy as np





def linbinning(values, xmin=-1, xmax=-1, n=10, silent=False):
	'''
	Returns a LINEAR binning of the data, `values`.
	Launch as:

	logbinning(values, xmin=-1, xmax=-1, n=10)

	Output: ibin, binvalue, histo[ibin], histo[ibin]/len(values)

	'''
	
	#Parameters
	if(xmin==-1):
		xmin=values.min()
	if(xmax==-1):
		xmax=values.max()
	assert(xmin>0 and xmax>xmin)
		
	delta=np.double(xmax-xmin)/(n-1)
	histo=np.zeros(n)

	for val in values:
		ibin=int((val-xmin)/delta)
		histo[ibin]+=1

	print("Linear binning")
	output=np.ndarray((n,4))
	for ibin in range(n):
		output[ibin]=[ibin,xmin+ibin*delta, histo[ibin], histo[ibin]/len(values)]
		if not silent:
			print("LIN",  ibin,xmin+ibin*delta, histo[ibin], histo[ibin]/len(values))

	return output


def logbinning(values, xmin=-1, xmax=-1, n=10, silent=False):
	'''
	Returns a LOGARITHMIC binning of the data, `values`.
	Launch as:

	logbinning(values, xmin=-1, xmax=-1, n=10)

	Make sure there are no zeros in `values`.

	Output: ibin, binvalue, histo[ibin], histo[ibin]/len(values)

	'''

	#Parameters
	if(xmin==-1):
		xmin=values.min()
	if(xmax==-1):
		xmax=values.max()
	assert(xmin>0 and xmax>xmin)
		
	#Grandezze derivate
	ymin=np.log(xmin)
	ymax=np.log(xmax+0.01*(xmax-xmin/n))
	delta=np.double(ymax-ymin)/n
	histo=np.zeros(n)

	for val in values:
		if val<xmin: continue
		yi=np.log(val)
		ibin=int((yi-ymin)/delta)
		if ibin<n:
			histo[ibin]+=1

	print("Logarithmic binning")
	output=np.ndarray((n,4))
	for ibin in range(len(histo)):
		output[ibin]=[ibin, np.exp(ymin+ibin*delta), histo[ibin], histo[ibin]/len(values)]
		if not silent:
			print("LOG",  ibin, np.exp(ymin+ibin*delta), histo[ibin], histo[ibin]/len(values))
		
	return output


if __name__=='__main__':
	nvalues=500
	expList=np.random.exponential(1,nvalues)
	linbinning(expList)
	logbinning(expList)
