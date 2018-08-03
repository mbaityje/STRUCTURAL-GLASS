import numpy as np
import matplotlib.pyplot as plt
from beylkinold import Beylkin
import sys

b = Beylkin()

f = np.load('CFF.npy')

t = np.arange(len(f))

b.driver_load(t, f, int(sys.argv[1]))

plt.plot(f)

plt.plot(b.correction())

plt.show()
