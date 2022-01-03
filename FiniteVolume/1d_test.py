import numpy as np
from numpy.linalg import norm
import sys
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from enum import Enum
class Grid:
    def __init__(self, lx, ly, lz, ncv):
        self._ncv = ncv
        dx = lx/float(ncv)
        self._xf = np.array([i*dx for i in range(ncv+1)])
        self._xP = np.array([self._xf[0]] +
                             [0.5*(self._xf[i]+self._xf[i+1]) for i in range(ncv)]
        +
                             [self._xf[-1]])
        self._Af = ly*lz*np.ones(ncv+1)
        self._Ao = (2.0*dx*ly + 2.0*dx*lz)*np.ones(ncv)
        self._vol = dx*ly*lz*np.ones(ncv)
    @property
    def ncv(self):
        return self._ncv
    @property
    def xf(self):
        return self._xf
    @property
    def xP(self):
        return self._xP
    @property
    def dx_WP(self):
        return self.xP[1:-1]-self.xP[0:-2]
    @property
    def dx_PE(self):
        return self.xP[2:]-self.xP[1:-1]
    @property
    def Af(self):
        return self._Af
    @property
    def Aw(self):
        return self._Af[0:-1]
    @property
    def Ae(self):
        return self._Af[1:]
    @property
    def Ao(self):
        return self._Ao
    @property
    def vol(self):
        return self._vol

class ScalarCoeffs:
    def __init__(self, ncv):
        self._ncv = ncv
        self._aP = np.zeros(ncv)
        self._aW = np.zeros(ncv)
        self._aE = np.zeros(ncv)
        self._rP = np.zeros(ncv)
    def zero(self):
        self._aP.fill(0.0)
        self._aW.fill(0.0)
        self._aE.fill(0.0)
        self._rP.fill(0.0)
    def accumulate_aP(self, aP):
        self._aP += aP
    def accumulate_aW(self, aW):
        self._aW += aW
    def accumulate_aE(self, aE):
        self._aE += aE
    def accumulate_rP(self, rP):
        self._rP += rP
    @property
    def ncv(self):
        return self._ncv
    @property
    def aP(self):
        return self._aP
    @property
    def aW(self):
        return self._aW
    @property
    def aE(self):
        return self._aE
    @property
    def rP(self):
        return self._rP

lx = 1.0
ly = 0.1
lz = 0.1
ncv = 10
grid = Grid(lx, ly, lz, ncv)
maxIter = 100
converged = 1e-6
k = 0.1
coeffs = ScalarCoeffs(grid.ncv)
print("Hello")
