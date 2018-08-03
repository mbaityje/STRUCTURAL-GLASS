#!/usr/bin/env python

"""
test_beylkin.py
Usage: python example.py filename.npy
Author: Ian S. Dunn, Columbia University
"""

import numpy as np
import matplotlib.pyplot as plt
from beylkin import Beylkin
import sys
filename='CFF.npy'

b = Beylkin(decaying=True)

f = np.load(filename)

t = np.arange(len(f))

b.driver_load(t, f, int(sys.argv[1]))

plt.plot(f,label=filename)

plt.plot(b.correction(),'o',markersize=3.,label='correction')

tt = np.arange(1 * len(f)) * 1 #a finer grid
#plt.plot(tt, b.prony_function(tt),label='prony_function')
plt.legend()
plt.show()
