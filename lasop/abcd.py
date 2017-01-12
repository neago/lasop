# -*- coding: utf-8 -*-
"""
Created on Mon Jul 05 18:55:55 2010
Last changed on 

@author: Jonas Neergaard-Nielsen

Tools for calculation of Gaussian beam propagation using ABCD matrix formalism.
All lengths are in units of mm.
"""

import numpy as np
from numpy import pi, conj

# set wavelength to 860 nm
lam = 0.000860

# ====================
# q-parameter formulas
#

def wR2q(w, R, n=1):
    """
    q = wR2q(w, R, n=1)
    --------------
    Get the q-parameter from a given spot size and radius of curvature.
    n is the medium's refractive index.
    """
    return 1/(1/R - 1j * lam/n / (pi * w**2))
    
    
def w02q(w0, n=1):
    """
    q = w02q(w0, n=1)
    ------------
    Get the q-parameter at a waist point from the waist size.
    n is the medium's refractive index.
    """
    return 1j * pi * w0**2 / (lam/n)


def q2w(q, n=1):
    """
    w = q2w(q, n=1)
    ----------
    Get the spot size from a given q-parameter.
    n is the medium's refractive index.
    """
    return np.sqrt(-lam/n / (pi * np.imag(1 / q)))


def q2R(q):
    """
    w = q2R(q, n=1)
    ----------
    Get the beam radius of curvature from a given q-parameter.
    n is the medium's refractive index.
    """
    return 1/ np.real(1 / q)
    

def q2w0(q, n=1):
    """
    w0 = q2w0(q, n=1)
    ------------
    Get the waist size from a given q-parameter.
    n is the medium's refractive index.
    """
    return np.sqrt(np.imag(q) * lam/n / pi)
    
    
def q2div(q, n=1):
    """
    div = q2div(q, n=1)
    --------------
    Get the far-field beam divergence for a given q-parameter.
    n is the medium's refractive index.
    """
    return lam/n / (pi * q2w0(q))
    
    
def qABCD(q, M):
    """
    q1 = qABCD(q0, M)
    -----------------
    Transform the q-parameter according to the ABCD matrix M.
    """
    M = np.array(M)
    return (M[0, 0] * q + M[0, 1]) / (M[1, 0] * q + M[1, 1])
    
    
def qreverse(q):
    """
    q1 = qreverse(q)
    ----------------
    q-parameter transformation when changing propagation direction.
    """
    return -conj(q)
    
    
def qpropagate(zini, qini, elements, z):
    """
    qout = qpropagate(zini, qini, elements, z)
    ------------------------------------------
    Propagate the q-parameter through an optical system.
    zini, qini : location and value of a known q-parameter of the beam
                 (qini must be given for forward propagation of the beam)
    elements   : list of [z-location, ABCD matrix] descriptions of the 
                 optical elements
    z          : location to calculate output q-parameter (if z < zini, the
                 output q-parameter will still be for forward propagation)
    """
    elements = elements[:]
    elements.sort()
    zt = zini
    qt = qini
    
    if z >= zini:
        elements.reverse()
        while elements:
            el = elements.pop()
            if zt <= el[0] <= z:
                qt += el[0] - zt
                qt = qABCD(qt, el[1])
                zt = el[0]
        qt += z - zt
    else:
        qt = qreverse(qt)
        while elements:
            el = elements.pop()
            if z <= el[0] <= zt:
                qt += zt - el[0]
                qt = qABCD(qt, el[1])
                zt = el[0]
        qt += zt - z
        qt = qreverse(qt)
    
    return qt
        
    
# =============
# ABCD matrices
#

def Mprop(d):
    """
    M = Mprop(d)
    ------------
    ABCD matrix for free space propagation of distance d.
    """
    return np.matrix([[1, d], [0, 1]])


def Minterface(n0, n1, R=np.inf):
    """
    M = Minterface(n0, n1, R='inf')
    ----------------------
    ABCD matrix for the refraction at an interface (with radius of curvature R) 
    from a medium with refractive index n0 to a medium with refractive index n1.
    If no R is given, R=infinite i.e. flat surface is assumed.
    R>0 means convex interface.
    """
    return np.matrix([[1, 0], [(n0-n1)/(R*n1), n0/n1]])


def Mlens(f):
    """
    M = Mlens(f)
    ------------
    ABCD matrix for a thin lens of focal length f.
    """
    return np.matrix([[1, 0], [-1/f, 1]])
    
    
def Mmirror(R):
    """
    M = Mmirror(R)
    --------------
    ABCD matrix for a curved mirror with radius of curvature R.
    Concave mirrors have R<0, convex have R>0.
    """
    return np.matrix([[1, 0], [2/R, 1]])


# ########################
# default initializations

