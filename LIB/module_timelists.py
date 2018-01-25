#!/usr/bin/python
#
# Module for time lists
# 
# 
from __future__ import print_function
import numpy as np


def ListaLogaritmica(x0,xn,n,ints=False,addzero=False):
    """
    Creates a list of logarithmically spaced times (or whatever you want them to be).
    output = ListaLogaritmica(x0,xn,n,ints=False,addzero=False)
    x0: initial time
    xn: final time
    n:  number of times
    ints: if True, the output takes only integer values. Every integer only appears once. Default is False.
    addzero: prepends the value 0 to the list
    """
    assert(xn>x0)
    assert(x0>0)
    n=np.int64(n)
    y0=np.log(x0)
    yn=np.log(xn)
    delta=np.float64(yn-y0)/(n-1)
    listax=np.exp([y0+i*delta for i in range(n)])
    if ints:
        listax=np.unique(np.round(listax,0).astype(int))
    if addzero:
        listax=np.insert(listax,0,0)
    return listax

