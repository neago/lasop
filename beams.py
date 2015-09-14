from pylab import *
from optics import abcd


class Beam():
    def __init__(self, w0=1, z0=0):
        self.w0 = w0
        self.z0 = z0
        self.q0 = abcd.w02q(w0)

    def width(self, z):
        return abcd.q2w(self.q0 + z)
        
    def plot_profile(self, zmin=0, zmax=100):
        z = linspace(zmin, zmax, 200)
        plot(z, self.width(z))


class Element:


class System():
    