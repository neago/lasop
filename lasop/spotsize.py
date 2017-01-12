# -*- coding: utf-8 -*-
"""
Created on Mon Mar 08 18:58:29 2010

@author: Jonas S. Neergaard-Nielsen

Script to load ASCII waveform trace from the LeCroy and fit it to a double
error function to estimate beam spot size.
"""

import os
import pylab
import scipy.optimize
from scipy.special import erf


def select_file(datadir, manual_choice=False):
    """
    Select single file among all files in datadir.
    If manual_choice = False, last file in dir will be used (sorted by 
    filename).
    """
    
    filelist = [os.path.join(datadir, fn) for fn in os.listdir(datadir)]
    filelist.sort()
    
    if manual_choice:
#        for i in range(len(filelist)):
#            print '[', str(i), '] ', os.path.basename(filelist[i])
        indices = '[0 - ' + str(len(filelist)-1) + ']'
        reply = input('Please choose a file index ' + indices
                          + ' (q to skip): ')
        if reply == 'q':
            return None
        else:
            index_input = int(reply)
            chosen = filelist[index_input]
    else:
        chosen = filelist[-1]
        
    return chosen


def read_lcascii(filename):
    """
    Read a single waveform comma-separated ASCII file from LeCroy 
    WaveRunner and return a list of (time, ampl) pairs.
    """
    
    x = []
    y = []
    flag = 0
    fin = open(filename)
    
    for line in fin:
        if flag:
            co = [float(val) for val in line.split(',')]
            x.append(co[0])
            y.append(co[1])
        else:
            if line == 'Time,Ampl\n': 
                flag = 1
            
    return x, y


def fit(x, y, p0=None):
    """
    Fit the (x, y) data to the erf-based function for a ruler passing
    through the beam.
    Returns (fitfunc(p, x), fitparam).
    """
    
    if p0 is None:    # initial parameter estimates, if not manually provided
        p0 = [min(y), max(y)-min(y), 0., (x[-1]-x[0])/50, (x[-1]-x[0])/2]
    
    fitfunc = lambda p, x: p[0] + p[1]*(1 - erf((x-p[2])/p[3])/2
                                        + erf((x-p[2]-p[4])/p[3])/2)
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    
    fitparam, success = scipy.optimize.leastsq(errfunc, p0[:], args=(x, y))
    
    if success:
        fitted_w = RULER_WIDTH * fitparam[3] / fitparam[4]
        return fitfunc, fitparam, fitted_w
    else:
        print("Fit error.")
        return None
        

def plot_fit(x, y, fitfunc, fitparam):
    """Plots data and fitted function."""
    
    pylab.figure()
    pylab.plot(x, y, 'b.')
    pylab.plot(x, fitfunc(fitparam, x), 'r-')
    
    fitted_w = RULER_WIDTH * fitparam[3] / fitparam[4]
    pylab.title('{0:.3g}'.format(fitted_w), fontsize = 42)


if __name__=='__main__':

    RULER_WIDTH = 12.7
    datadir = input('Directory? ')
        
    while True:
        fn = select_file(datadir, True)
        if fn:
            x, y = read_lcascii(fn)
            input_params = input('Manually input parameters? ' +
                                     '(Enter to use defaults): ')
            if input_params:
                test_params = input('Input comma-separated values:\n' +
                                        '0-level, ampl., left edge offset, w, ' +
                                        'ruler width\n')
                test_params = [float(p) for p in test_params.split(',')]
            else:
                test_params = None
                
            fitfunc, fitparam, fitted_w = fit(x, y, test_params)
            plot_fit(x, y, fitfunc, fitparam)
            print((os.path.basename(fn) + ': ' + '{0:.4g}'.format(fitted_w) + ' mm'))
        
        stopnow = input('Press Enter for new file, q to quit: ')
        if stopnow == 'q':
            break



