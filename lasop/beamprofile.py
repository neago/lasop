# -*- coding: utf-8 -*-
"""
Created on Fri Jul 02 11:03:29 2010

@author: exp
"""

import os
#import numpy as np
#import scipy as sp
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import scipy.optimize
#from scipy.special import erf

from numpy import pi, sqrt
from numpy import arange
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

lam = .000860     # wavelength in mm

def sslineparse(line):
    if line[0] == '#':
        return
    
    line = line.split('\t')
    point = [eval(x) for x in line]
    return point
    

def ssfileparse(file):
    lines = [ln for ln in open(file)]
    data = [sslineparse(ln) for ln in lines if sslineparse(ln)]
    x = [p[0] for p in data]
    y = [p[1] for p in data]
    return x, y
    
def fit(x, y, p0=None):
    
    if p0 is None:
        p0 = [.1, 0]
        
    fitfunc = lambda p, z: p[0] * sqrt(1 + ((z - p[1]) / (pi*p[0]**2 / lam))**2)
    errfunc = lambda p, z, y: fitfunc(p, z) - y
    
    fitparam, success = leastsq(errfunc, p0[:], args=(x, y))
    
    if not success:
        print 'Fit error'
        return
    return fitfunc, fitparam


# script code
if __name__ == "__main__":
    datadir = "d:/data/categorized/spotsize"
    filelist = [os.path.join(datadir, fn) for fn in os.listdir(datadir)]
    
    for file in filelist:
        print filelist.index(file), '\t', file
    fi = input('Which file? ')
        
    x, y = ssfileparse(filelist[fi])
    func, param = fit(x, y, [.1,100])
    
    plt.plot(x, y, 'o')
    pxmin = min(0, min(x))
    pxmax = max(0, max(x))
    pxrange = pxmax - pxmin
    ptpoints = 50
    px = arange(pxmin - pxrange/10., pxmax + pxrange/10., 1.2*pxrange/ptpoints)
    plt.plot(px, func(param, px), 'r-')
    print "w_0 = %.4f mm;  z_0 = %.1f mm" % tuple(param)