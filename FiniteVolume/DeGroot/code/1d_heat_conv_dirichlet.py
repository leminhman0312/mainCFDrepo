# Solving the 1D Diffusion Equation
import numpy as np
from enum import Enum
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from numpy.linalg import norm
import sys
import matplotlib.pyplot as plt


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
        self._xP = np.array([self._xf[0]] +
                             [0.5*(self._xf[i]+self._xf[i+1]) for i in range(ncv)] +
                             [self._xf[-1]])
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

    def accumulate_rP(self,rP):
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
        self._gradient = gradient
        self._loc = loc

    def value(self):
        # return the boundary condition value
        if self._loc is BoundaryLocation.WEST:
            return self._phi[1]-self._gradient*self._grid.dx_WP[0] # T(0) = T(1) - gb*dx
        elif self._loc is BoundaryLocation.EAST:
            return self._phi[-2]-self._gradient*self._grid.dx_PE[-1] #T(-1) = T(-2) -gb*dx
        else:
            raise ValueError("Unknown boundary location")


    def coeff(self):
        # return the linearization coefficient
        return 1

    def apply(self):
        # applies the boundary condition in the referenced field variable array
        if self._loc is BoundaryLocation.WEST:
            self._phi[0] = self._phi[1] - self._gradient*self._grid.dx_WP[0]
        elif self._loc is BoundaryLocation.EAST:
            self._phi[-1] = self._phi[-2] + self._gradient*self._grid.dx_PE[-1]
        else:
            raise ValueError("Unknown boundary location")

class RobinBC:
    def __init__(self,phi,grid,h,k,tinfty, loc):
        # # constructor
        # phi = field variable array
        # grid = grid
        # h = convective coefficient
        # k = conductive coefficient
        # tinfty = freestream temperature
        # loc = boundary location

        self._phi = phi
        self._grid = grid
        self._h = h
        self._k = k
        self._tinfty = tinfty
        self._loc = loc

    def value(self):
        # return the boundary condition value
        if self._loc is BoundaryLocation.WEST:
            return (self._phi[1] + self._grid.dx_WP[0]*(self._h/self._k)*(self._tinfty))/(1 + self._grid.dx_WP[0]*(self._h/self._k))
        elif self._loc is BoundaryLocation.EAST:
            return (self._phi[-2] - self._grid.dx_PE[-1]*(self._h/self._k)*(self._tinfty))/(1 - self._grid.dx_PE[-1]*(self._h/self._k))
        else:
            raise ValueError("Unknown boundary location")


    def coeff(self):
        # return the linearization coefficient
        return 1

    def apply(self):
        # applies the boundary condition in the referenced field variable array
        if self._loc is BoundaryLocation.WEST:
            self._phi[0]  =  (self._phi[1] + self._grid.dx_WP[0]*(self._h/self._k)*(self._tinfty))/(1 + self._grid.dx_WP[0]*(self._h/self._k))
        elif self._loc is BoundaryLocation.EAST:
            self._phi[-1] = (self._phi[-2] - self._grid.dx_PE[-1]*(self._h/self._k)*(self._tinfty))/(1 - self._grid.dx_PE[-1]*(self._h/self._k))
        else:
            raise ValueError("Unknown boundary location")



    

# class DiffusionModel defining a defusion model
class DiffusionModel:
    def __init__(self,grid,phi,gamma,west_bc,east_bc):
        # constructor
        self._grid = grid
        self._phi = phi
        self._gamma = gamma     # gamma here is "k" in the heat model
        self._west_bc = west_bc
        self._east_bc = east_bc



    # add diffusion terms to coefficient arrays
    def add(self,coeffs):

        # calculate west/east diffusion flux for each face
        flux_w = -self._gamma*self._grid.Aw*(self._phi[1:-1]-self._phi[0:-2])/self._grid.dx_WP # Fw = k*(Tp-Tw)*Aw/ dx_WP
        flux_e = -self._gamma*self._grid.Ae*(self._phi[2:]-self._phi[1:-1])/self._grid.dx_PE # Fe = k*(Te-Tp)*Ae/dx_PE

        # calculate the linearized coefficient [aW, aP, aE]
        coeffW = -self._gamma*self._grid.Aw/self._grid.dx_WP
        coeffE = -self._gamma*self._grid.Ae/self._grid.dx_PE
        coeffP = -coeffW - coeffE # should be +?

        # modified the linearized coefficient on the boundaries
        coeffP[0] += coeffW[0]*self._west_bc.coeff()
        coeffP[-1] += coeffE[-1]*self._east_bc.coeff()

        # zero the coefficients that are not used (on the boundary)
        coeffW[0] = 0.0
        coeffE[-1] = 0.0

        # calculate the net flux from each cell (out -in)
        flux = flux_e - flux_w

        # add to coefficient arrays
        coeffs.accumulate_aP(coeffP)
        coeffs.accumulate_aW(coeffW)
        coeffs.accumulate_aE(coeffE)
        coeffs.accumulate_rP(flux)

        # return the coefficient arrays
        return coeffs


# function to return a sparse matrix representation of a set of scalar coefficients
def get_sparse_matrix(coeffs):
    ncv = coeffs.ncv
    data = np.zeros(3*ncv-2)
    rows = np.zeros(3*ncv-2, dtype = int)
    cols = np.zeros(3*ncv-2, dtype = int)
    data[0] = coeffs.aP[0]
    rows[0] = 0
    cols[0] = 0

    if ncv > 1:
        data[1] = coeffs.aE[0]
        rows[1] = 0
        cols[1] = 1

    for i in range(ncv-2):
        data[3*i+2] = coeffs.aW[i+1]
        data[3*i+3] = coeffs.aP[i+1]
        data[3*i+4] = coeffs.aE[i+1]

        rows[3*i+2:3*i+5] = i+1

        cols[3*i+2] = i
        cols[3*i+3] = i+1
        cols[3*i+4] = i+2

    if ncv > 1:
        data[3*ncv-4] = coeffs.aW[-1]
        data[3*ncv-3] = coeffs.aP[-1]
        
        rows[3*ncv-4:3*ncv-2] = ncv-1

        cols[3*ncv-4] = ncv-2
        cols[3*ncv-3] = ncv-1

    return csr_matrix((data, (rows,cols)))
                     
# solve the linear system and return the field variables
def solve(coeffs):
    # get the sparse matrix
    A = get_sparse_matrix(coeffs)
    # solve the linear system
    return spsolve(A, -coeffs.rP)
    



# class defining a surface convection model
class SurfaceConvectionModel:
    # constructor
    def __init__(self,grid,T,ho,To):
        self._grid = grid
        self._T = T
        self._ho = ho
        self._To = To

    # add surface convection terms to coefficient arrays
    def add(self,coeffs):
        # calculate the source term
        source = self._ho * self._grid.Ao*(self._T[1:-1]-self._To)

        # calculat linearization coefficients
        coeffP = self._ho * self._grid.Ao

        # add to coefficient arrays
        coeffs.accumulate_aP(coeffP)
        coeffs.accumulate_rP(source)

        return coeffs
        
        


# Solving the 1D steady conduction with Dirichet BC
# Initially, east = 300K, west = 300K
# define the grid
lx = 1.0
ly = 0.1
lz = 0.1
ncv = 10
mygrid = Grid(lx,ly,lz,ncv)

#set the max iterations and convergence criterion
maxIter = 100
converged=1e-6

# thermal properties
k = 100

# define convection parameters
ho = 25
To = 200

# coefficients
coeffs = ScalarCoeffs(mygrid.ncv)

# # initial condition
T0 = 300

# # initialize field variable arrays
T = T0*np.ones(mygrid.ncv+2)

# # # boundary condition
west_bc = DirichletBC(T, mygrid, 400, BoundaryLocation.WEST)
# east_bc = DirichletBC(T, mygrid, 300, BoundaryLocation.EAST)
east_bc = NeumannBC(T, mygrid, 0, BoundaryLocation.EAST)

# # # apply boundary conditions
west_bc.apply()
east_bc.apply()


# # # list to store the solution at each iteration
T_solns = [np.copy(T)]

# # # define the diffusion model
diffusion = DiffusionModel(mygrid, T, k, west_bc, east_bc)

# define the surface convection model
surfaceConvection = SurfaceConvectionModel(mygrid,T,ho,To)

avgResidList = []
iterList = []
# # # iterate until the solution is converged
for i in range(maxIter):

    # zero the coeff and add 
    coeffs.zero()
    coeffs = diffusion.add(coeffs)
    coeffs = surfaceConvection.add(coeffs)

    # compute residual and check for convergence
    maxResid = norm(coeffs.rP, np.inf)
    avgResid = np.mean(np.absolute(coeffs.rP))
    print("Iteration = {} ; Max Resid = {} ; Avg Resid = {}".format(i,maxResid, avgResid) )
    if maxResid < converged:
        break

    # solve the sparse matrix system
    dT = solve(coeffs)

    # update the solution and boundary conditions
    T[1:-1] += dT
    west_bc.apply()
    east_bc.apply()

    # store the solution
    T_solns.append(np.copy(T))
        
    avgResidList.append(avgResid)
    iterList.append(i)


# plotting
# i = 0
# for Ti in T_solns:
#     plt.plot(mygrid.xP, Ti, label = str(i))
#     i += 1
    
# plt.xlabel('X')
# plt.ylabel('T')
# plt.legend(loc='right')
# plt.show()

# print(iterList)

# plt.title('Residual')
# plt.plot(avgResidList,iterList)
# plt.show()
