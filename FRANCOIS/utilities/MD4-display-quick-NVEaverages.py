import numpy as np
import sys

filename = sys.argv[1]

things=np.loadtxt(filename)

averages = np.mean(things,0)
print averages
print 'avg. temperature:', averages[1]
print 'avg. pressure:', averages[2]
print 'avg. Ep:', averages[3]
print 'avg. Ec:', averages[4]
#print 'avg. momentum:', averages[5]
#print 'avg. volum:', averages[6]
#print 'avg. Natoms:', averages[7]


stds = np.std(things,0)
print stds
print 'std. dev. temperature:', stds[1]
print 'std. dev. pressure:', stds[2]
print 'std. dev. Ep:', stds[3]
print 'std. dev. Ec:', stds[4]

