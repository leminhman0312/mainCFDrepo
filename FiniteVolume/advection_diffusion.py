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
        flux_w = -self._gamma*self._grid.Aw*(self._phi[1:-1]-self._phi[0:-2])\
            /self._grid.dx_WP # Fw = k*(Tp-Tw)*Aw/ dx_WP
        flux_e = -self._gamma*self._grid.Ae*(self._phi[2:]-self._phi[1:-1])\
            /self._grid.dx_PE # Fe = k*(Te-Tp)*Ae/dx_PE

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

# defining a first order implicit transient model
class FirstOrderTransientModel:
    # constructor
    def __init__(self,grid,T,Told,rho,cp,dt):
        self._grid = grid
        self._T = T
        self._Told = Told
        self._rho = rho
        self._cp = cp
        self._dt = dt

    def add(self,coeffs):
        # calculate the transient term
        # rp  = rho cp Vp (TP - TPold)/ dt
        transient = self._rho*self._cp*self._grid.vol*(self._T[1:-1]-self._Told[1:-1])/self._dt
        
        # calculate the linearization coefficients
        coeff = self._rho*self._cp*self._grid.vol/self._dt # ap = rho*cp*vp/dt

        # add to coefficient arrays
        coeffs.accumulate_aP(coeff)
        coeffs.accumulate_rP(transient)

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
    
# class defining an upwind advection model
class UpwindAdvectionModel:

    # constructor
    def __init__(self,grid,phi,Uhe, rho, cp, west_bc, east_bc):
        self._grid = grid
        self._phi = phi
        self._Uhe = Uhe
        self._rho = rho
        self._cp = cp
        self._west_bc = west_bc
        self._east_bc = east_bc
        self._alphae = np.zeros(self._grid.ncv+1)
        self.phie = np.zeros(self._grid.ncv+1)


    # add diffusion terms to coefficient arrays
    def add(self,coeffs):

        # calculate the weighting factor
        for i in range(self._grid.ncv+1):
            if self._Uhe[i] >= 0:
                self._alphae[i] = 1
            else:
                self._alphae[i] = -1

        # calculate the east integration point (including both boundaries)
        self._phie= (1 + self._alphae)/2*self._phi[0:-1] + (1 - self._alphae)/2*self._phi[1:]

        # calculate the face mass fluxes
        mdote = self._rho*self._Uhe*self._grid.Af

        # calculate west/east advection flux
        flux_w = self._cp*mdote[:-1]*self._phie[:-1]
        flux_e = self._cp*mdote[1:]*self._phie[:-1]

        # calculate mass imbalance
        imbalance = -self._cp*mdote[1:]*self._phi[1:-1] + self._cp*mdote[:-1]*self._phi[1:-1]

        # calculate linearization coefficients
        coeffW = -self._cp*mdote[:-1]*(1 + self._alphae[:-1])/2
        coeffE = self._cp*mdote[1:]*(1 + self._alphae[1:])/2
        coeffP = -coeffW - coeffE

        # modify linearization coefficients on the boundaries
        coeffP[0] += coeffW[0] * self._west_bc.coeff()
        coeffP[-1] += coeffE[-1] * self._east_bc.coeff()

        # zero the bc that are not used
        coeffW[0] = 0.0
        coeffE[-1] = 0.0

        # calculate the net flux from each cell
        flux = flux_e -flux_w
        

        # add to coefficient arrays
        coeffs.accumulate_aE(coeffP)
        coeffs.accumulate_aW(coeffW)
        coeffs.accumulate_aE(coeffE)
        coeffs.accumulate_rP(flux)
        coeffs.accumulate_rP(imbalance)

        # return modified coeff arrays
        return coeffs

            
# Define the grid
lx = 1.0
ly = 0.1
lz = 0.1
ncv = 50
grid = Grid(lx,ly,lz,ncv)

# timestep information
nTime = 1
dt = 1e9
time = 0

# set iteration and convergence criterion
maxIter = 10
converged = 1e-6

# thermophysical properties
rho = 1000
cp = 4000
k = 0.5

# surface convection parameters
ho = 50
To = 200

# define coefficients
coeffs = ScalarCoeffs(grid.ncv)

# Initial conditions
T0 = 300
U0 = 0.01


# Initialize field variable arrays
T = T0*np.ones(grid.ncv+2)
Uhe = U0*np.ones(grid.ncv+1)

# define boundary conditions
west_bc = DirichletBC(T,grid,400,BoundaryLocation.WEST)
east_bc = DirichletBC(T,grid,0,BoundaryLocation.EAST)

# apply boundary conditions
west_bc.apply()
east_bc.apply()

# define the transient model
Told = np.copy(T)
transient = FirstOrderTransientModel(grid,T,Told,rho,cp,dt)

# define the diffusion model
diffusion = DiffusionModel(grid,T,k,west_bc,east_bc)

# define the surface convection model
surfaceConvection = SurfaceConvectionModel(grid,T,ho,To)

# define the advection model
advection = UpwindAdvectionModel(grid,T,Uhe,rho,cp,west_bc,east_bc)


# loop through all timesteps
for tStep in range(nTime):
    # update time
    time += dt

    # print timestep informatin
    print("Timestep = {}; Time = {}".format(tStep,time))

    # store the old temperature field
    Told[:] = T[:]

    # iterate until solution is converged
    for i in range(maxIter):
        # zero out the coefficents and add
        coeffs.zero()
        coeffs = diffusion.add(coeffs)
        coeffs = surfaceConvection.add(coeffs)
        coeffs = advection.add(coeffs)
        coeffs = transient.add(coeffs)
        
    
        # compute residual and check for convergence
        maxResid = norm(coeffs.rP, np.inf)
        avgResid = np.mean(np.absolute(coeffs.rP))
        print("Iteration = {} ; MaxResid = {} ; AvgResid = {}".format(i,maxResid, avgResid) )
        
        if maxResid < converged:
            break

        # solve the sparse matrix system
        dT = solve(coeffs)

        # update the solution and boundary conditions
        T[1:-1] += dT
        west_bc.apply()
        east_bc.apply()


plt.plot(grid.xP, T)
plt.xlabel("X")
plt.ylabel("T")
plt.show()
