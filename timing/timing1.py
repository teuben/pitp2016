#! /usr/bin/env python
#
#
#  demo:   using "run -i timing1.py"

import numpy as np


def make1(n, mode=0):
    if mode==0:
        return np.zeros(n, dtype=np.float64)         # zeros
    elif mode==1:   
        return np.arange(float(n))                   # ordinals
    elif mode==2:
        return np.random.uniform(0.0,1.0,int(n))     # random, uniform
    elif mode==3:
        return np.random.normal(0.0,1.0,int(n))      # random, normal


def make2(n, mode=0):
    if mode==0:
        a=[]
        for i in range(n):
            a.append(0)
        return a
    elif mode==1:
        a=[]
        for i in range(n):
            a.append(i)
        return a
    elif mode==2:
        return list(range(n))       # for python3 , otherwise this is a constant independant of n !!

def delta1(a):
    n = len(a)
    d = np.zeros(n-1)
    for i in range(1,n):
        d[i-1] = a[i]-a[i-1]
    return d


def delta2(a):
    return a[1:] - a[:-1]

    
