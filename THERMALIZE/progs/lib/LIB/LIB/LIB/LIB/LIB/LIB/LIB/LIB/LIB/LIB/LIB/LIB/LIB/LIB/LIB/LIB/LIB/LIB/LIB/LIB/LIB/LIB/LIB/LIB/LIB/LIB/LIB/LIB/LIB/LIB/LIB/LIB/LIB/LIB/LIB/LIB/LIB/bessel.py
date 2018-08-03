
import numpy as np
from scipy.special import jv

def bessel(x):

    b = 100. * np.pi
    
    return jv(0, b * x)
