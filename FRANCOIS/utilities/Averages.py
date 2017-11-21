#
# Takes several text files as input, and calculates mean and standard error
# of each. All files must have the same format.
#
# INPUT:
# a b c ....
#
# OUTPUT
# a_mean b_mean c_mean ... a_err b_err c_err
#
import sys
import numpy as np
from pprint import pprint

if len(sys.argv)==1:
    print("Launch as:")
    print("python Averages.py file1 file2 file3 ...")
    print("OUTPUT: first, all the averages of each column, then, all the errorbars")
    sys.exit()

data = []
my_files=sys.argv[1:]
print(my_files)
for filename in my_files:
    a = np.array(np.loadtxt(filename, dtype=np.float64, delimiter = ' ', skiprows = 0))
    data.append(a)
data = np.array(data)

media=np.mean(data,dtype=np.float64,axis=0)
devst=np.std(data,dtype=np.float64,axis=0)/np.float64((len(data)-1))
all=np.concatenate((media,devst),axis=1)

np.savetxt(sys.stdout,all)




