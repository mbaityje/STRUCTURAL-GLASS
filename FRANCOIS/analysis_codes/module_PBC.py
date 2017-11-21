#module_PBC
from __future__ import print_function
import numpy as np
from numba import jit
from numba import float32, int32, int64, float64


## assume isotropy : same PBC in all directions ! ##
@jit(nopython=True, nogil=True)
def applyPBC_vector_func(displacements, Lx):
    displacements[(displacements >  Lx/2.0)] -= Lx
    displacements[(displacements <- Lx/2.0)] += Lx
    return displacements    # we want the coreectly signed displacement here. 

@jit(nopython=True, nogil=True)
def applyPBC_pair_func(x1, x2, Lx):
    dist = x1-x2
    if   dist > Lx/2.0:
        dist -= Lx
    elif dist < -Lx/2.0:
        dist += Lx
    return dist
