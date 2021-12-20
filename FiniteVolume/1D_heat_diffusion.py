# Solving the 1D Diffusion Equation
import numpy as np
from enum import Enum


# class Grid: defining a 1D cartesian grid
class Grid:

    def __init__(self,lx,ly,lz,ncv):
        # # constructor
        # lx = total length of domain in x direction [m]
        # ly = total length of domain in y direction [m]
        # lz = total length of domain in z direction [m]

        # store the number of control volumes
        self._ncv = ncv

        # calculate the control volume length
        dx = lx/float(ncv)

        # calculate the face locations
        # generate array of control volume element with length dx each
        self._xf = np.array([i*dx for i in range(ncv+1)])

        # calculate the cell centroid locations
        # Note: figure out why add a xf[-1] last item 
        self._xP = np.array([self._xf[0] +
                             [0.5*(self._xf[i]+self._xf[i+1]) for i in range(ncv)] +
                             [self._xf[-1]]])
        # calculate face areas
        self._Af = ly*lz*np.ones(ncv+1)

        # calculate the outer surface area of each cell
        self._Ao = (2.0*dx*ly + 2.0*dx*lz)*np.ones(ncv)

        # calculate cell volumes
        self._vol = dx*ly*lz*np.ones(ncv)

        @property
        def ncv(self):
            # number of control volumes in domain
            return self._ncv

        @property
        def xf(self):
            # face location array
            return self._xf

        @property
        def xP(self):
            # cell centroid array
            return self._xP

        @property
        def dx_WP(self):
            return self.xP[1:-1]-self.xP[0:-2]

        @property
        def dx_PE(self):
            return self.xP[2:]-self.xP[1:-1]

        @property
        def Af(self):
            # Face area array
            return self._Af

        @property
        def Aw(self):
            # West face area array
            return self._Af[0:-1]

        @property
        def Ae(self):
            # East face area array
            return self._Af[1:]

        @property
        def Ao(self):
            # Outer face area array
            return self._Ao

        @property
        def vol(self):
            # Cell volume array
            return self._vol



# class ScalarCoeffs defining coeficients for the linear system
class ScalarCoeffs:
    def __init__(self,ncv):
        # # constructor:
        # ncv = number of control volume in domain

        self._ncv = ncv
        self._aP = np.zeros(ncv)
        self._aW = np.zeros(ncv)
        self._aE = np.zeros(ncv)
        self._rP = np.zeros(ncv)

    def zero(self):
        # zero out the coefficient arrays
        self._aP.fill(0.0)
        self._aW.fill(0.0)
        self._aE.fill(0.0)
        self._rP.fill(0.0)

    def accumulate_aP(self,aP):
        # accumulate values onto aP
        self._aP += aP

    def accumulate_aW(self,aW):
        # accumulate values onto aW
        self._aW += aW
    
    def accumulate_aE(self,aE):
        # accumulate values onto aE
        self._aE += aE

    def accumulate_aP(self,rP):
        # accumulate values onto rP
        self._rP += rP
    
    @property
    def ncv(self):
        # number of control volume in domain
        return self._ncv

    @property
    def aP(self):
        # cell coefficient
        return self._aP

    @property
    def aW(self):
        # West cell coefficient
        return self._aW
    
    @property
    def aE(self):
        # East cell coefficient
        return self._aE

    @property
    def rP(self):
        # cell coefficient
        return self._rP
    

# class BoundaryLocation defining boundary condition locations
class BoundaryLocation(Enum):
    WEST = 1
    EAST = 2
    

# Implementing Boundary Conditions
class DirichletBC:
    def __init__(self,phi,grid,value,loc):
        # # constructor
        # phi = field variable array
        # grid = grid
        # value = boundary value
        # loc = boundary location

        self._phi = phi
        self._grid = grid
        self._value = value
        self._loc = loc

    def value(self):
        # return the boundary condition value
        return self._value

    def coeff(self):
        # return the linearization coefficient
        return 0

    def apply(self):
        # applies the boundary condition in the referenced field variable array
        if self._loc is BoundaryLocation.WEST:
            self._phi[0] = self._value
        elif self._loc is BoundaryLocation.EAST:
            self._phi[-1] = self._value
        else:
            raise ValueError("Unknown boundary location")


class NeumannBC:
    def __init__(self,phi,grid,gradient,loc):
        # # constructor
        # phi = field variable array
        # grid = grid
        # gradient = gradient at cell adjacent to boundary
        # loc = boundary location

        self._phi = phi
        self._grid = grid
        self._value = value
        self._loc = loc

    def value(self):
        # return the boundary condition value
        if self._loc is BoundaryLocation.WEST:
            return self._phi[1]-self._gradient*self._grid.dx_WP[0]
        elif self._loc is BoundaryLocation.EAST:
            return self._phi[-2]-self._gradient*self._grid.dx_WP[-1]
        else:
            raise ValueError("Unknown boundary location")


    def coeff(self):
        # return the linearization coefficient
        return 0

    def apply(self):
        # applies the boundary condition in the referenced field variable array
        if self._loc is BoundaryLocation.WEST:
            self._phi[0] = self._phi[1] - self._gradient*self._grid.dx_WP[0]
        elif self._loc is BoundaryLocation.EAST:
            self._phi[-1] = self._phi[-2] + self._gradient*self._grid.dx_WP[-1]
        else:
            raise ValueError("Unknown boundary location")
