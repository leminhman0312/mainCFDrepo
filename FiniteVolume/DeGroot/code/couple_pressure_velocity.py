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
    def __init__(self,grid,phi,phiold,rho,const,dt):
        self._grid = grid
        self._phi = phi
        self._phiold = phiold
        self._rho = rho
        self._const = const
        self._dt = dt

    def add(self,coeffs):
        # calculate the transient term
        # rp  = rho cp Vp (TP - TPold)/ dt
        transient = self._rho*self._const*self._grid.vol*(self._phi[1:-1]-self._phiold[1:-1])/self._dt
        
        # calculate the linearization coefficients
        coeff = self._rho*self._const*self._grid.vol/self._dt # ap = rho*cp*vp/dt

        # add to coefficient arrays
        coeffs.accumulate_aP(coeff)
        coeffs.accumulate_rP(transient)

        return coeffs


def get_sparse_matrix(PP_coeffs, PU_coeffs, UP_coeffs, UU_coeffs):
    """Function to return a sparse matrix representation of a set of scalar coefficients"""

    # Get number of control volumes (check that all are consistent)
    ncv = PP_coeffs.ncv
    if ncv is not PU_coeffs.ncv or ncv is not UP_coeffs.ncv or ncv is not UU_coeffs.ncv:
        raise ValueError("Not all coefficient arrays have the same dimension")

    # Set up the data and indexing arrays
    nvar = 2  # P, U
    ncoeff = 3  # aW, aP, and aE
    ndata = nvar*nvar*ncoeff*ncv - 2*nvar*nvar # Need to subtract boundary coefficients
    data = np.zeros(ndata)
    rows = np.zeros(ndata, dtype=int)
    cols = np.zeros(ndata, dtype=int)

    # Set up the first cell
    data[0] = PP_coeffs.aP[0]
    cols[0] = 0
    data[1] = PU_coeffs.aP[0]
    cols[1] = 1
    data[2] = PP_coeffs.aE[0]
    cols[2] = 2
    data[3] = PU_coeffs.aE[0]
    cols[3] = 3
    data[4] = UP_coeffs.aP[0]
    cols[4] = 0
    data[5] = UU_coeffs.aP[0]
    cols[5] = 1
    data[6] = UP_coeffs.aE[0]
    cols[6] = 2
    data[7] = UU_coeffs.aE[0]
    cols[7] = 3

    rows[0:4] = 0
    rows[4:8] = 1

    # Set up the interior cells
    for i in range(1, ncv-1):
        start = nvar*nvar*ncoeff*(i-1) + 8
        data[start+0] = PP_coeffs.aW[i]
        cols[start+0] = nvar*i - 2
        data[start+1] = PU_coeffs.aW[i]
        cols[start+1] = nvar*i - 1
        data[start+2] = PP_coeffs.aP[i]
        cols[start+2] = nvar*i
        data[start+3] = PU_coeffs.aP[i]
        cols[start+3] = nvar*i + 1
        data[start+4] = PP_coeffs.aE[i]
        cols[start+4] = nvar*i + 2
        data[start+5] = PU_coeffs.aE[i]
        cols[start+5] = nvar*i + 3
        data[start+6] = UP_coeffs.aW[i]
        cols[start+6] = nvar*i - 2
        data[start+7] = UU_coeffs.aW[i]
        cols[start+7] = nvar*i - 1
        data[start+8] = UP_coeffs.aP[i]
        cols[start+8] = nvar*i
        data[start+9] = UU_coeffs.aP[i]
        cols[start+9] = nvar*i + 1
        data[start+10] = UP_coeffs.aE[i]
        cols[start+10] = nvar*i + 2
        data[start+11] = UU_coeffs.aE[i]
        cols[start+11] = nvar*i + 3

        rows[start:start+6] = nvar*i
        rows[start+6:start+12] = nvar*i + 1

    # Set up the last cell
    i = ncv - 1
    start = nvar*nvar*ncoeff*(i-1) + 8
    data[start+0] = PP_coeffs.aW[i]
    cols[start+0] = nvar*i - 2
    data[start+1] = PU_coeffs.aW[i]
    cols[start+1] = nvar*i - 1
    data[start+2] = PP_coeffs.aP[i]
    cols[start+2] = nvar*i
    data[start+3] = PU_coeffs.aP[i]
    cols[start+3] = nvar*i + 1
    data[start+4] = UP_coeffs.aW[i]
    cols[start+4] = nvar*i - 2
    data[start+5] = UU_coeffs.aW[i]
    cols[start+5] = nvar*i - 1
    data[start+6] = UP_coeffs.aP[i]
    cols[start+6] = nvar*i
    data[start+7] = UU_coeffs.aP[i]
    cols[start+7] = nvar*i + 1

    rows[start:start+4] = nvar*i
    rows[start+4:start+8] = nvar*i + 1

    # Return the matrix
    return csr_matrix((data, (rows, cols)))

def solve(PP_coeffs, PU_coeffs, UP_coeffs, UU_coeffs):
    """Function to solve the linear system and return the correction fields"""

    # Get number of control volumes (check that all are consistent)
    ncv = PP_coeffs.ncv
    if ncv is not PU_coeffs.ncv or ncv is not UP_coeffs.ncv or ncv is not UU_coeffs.ncv:
        raise ValueError("Not all coefficient arrays have the same dimension")

    # Get the sparse matrix
    A = get_sparse_matrix(PP_coeffs, PU_coeffs, UP_coeffs, UU_coeffs)

    # Create and fill the right side vector
    b = np.zeros(2*ncv)
    for i in range(ncv):
        b[2*i] = - PP_coeffs.rP[i] - PU_coeffs.rP[i]
        b[2*i+1] = - UP_coeffs.rP[i] - UU_coeffs.rP[i]

    # Solve the linear system
    res = spsolve(A,b)

    # Extract the correction fields
    dP = np.zeros(ncv)
    dU = np.zeros(ncv)
    for i in range(ncv):
        dP[i] = res[2*i]
        dU[i] = res[2*i+1]

    return dP, dU




    
    
# class defining an upwind advection model
class UpwindAdvectionModel:

    # constructor
    def __init__(self,grid,phi,Uhe, rho, const, west_bc, east_bc):
        self._grid = grid
        self._phi = phi
        self._Uhe = Uhe
        self._rho = rho
        self._const = const
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
        flux_w = self._const*mdote[:-1]*self._phie[:-1]
        flux_e = self._const*mdote[1:]*self._phie[:-1]

        # calculate mass imbalance
        imbalance = -self._const*mdote[1:]*self._phi[1:-1] + self._const*mdote[:-1]*self._phi[1:-1]

        # calculate linearization coefficients
        coeffW = -self._const*mdote[:-1]*(1 + self._alphae[:-1])/2
        coeffE = self._const*mdote[1:]*(1 + self._alphae[1:])/2
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
        coeffs.accumulate_aP(coeffP)
        coeffs.accumulate_aW(coeffW)
        coeffs.accumulate_aE(coeffE)
        coeffs.accumulate_rP(flux)
        coeffs.accumulate_rP(imbalance)

        # return modified coeff arrays
        return coeffs

# class define advecting velocity model 
class AdvectingVelocityModel:
    # constructor
    def __init__(self,grid,dhat,Uhe,P,U,coeffs):
        self._grid = grid
        self._dhat = dhat
        self._Uhe  = Uhe
        self._P = P
        self._U = U
        self._coeffs = coeffs

    # update the advecting velocity array 
    def update(self):
        # calculate pressure gradient across faces
        gradPw = (self._P[1:-1]-self._P[0:-2])/self._grid.dx_WP
        gradPe = (self._P[2:]-self._P[1:-1])/self._grid.dx_PE

        # calculate cell pressure gradients
        gradP = 0.5*(gradPw + gradPe)

        # Calculating damping coefficient, dhat
        Ve = 0.5*(self._grid.vol[0:-1] + self._grid.vol[1:])
        ae = 0.5*(self._coeffs.aP[0:-1] + self._coeffs.aP[1:])
        self._dhat[1:-1] = Ve/ae

        # update advecting velocity
        self._Uhe[0] = self._U[0]
        self._Uhe[1:-1] = 0.5*(self._U[1:-2] + self._U[2:-1]) - self._dhat[1:-1]*(gradPe[:-1] - 0.5*(gradP[:-1] + gradP[1:]))
        self._Uhe[-1] = self._U[-1]

# class define pressure force model
class PressureForceModel:
    # constructor
    def __init__(self,grid,P,west_bc,east_bc):
        self._grid = grid
        self._P = P
        self._west_bc = west_bc
        self._east_bc = east_bc

    # function to add diffusion terms to coefficient arrays    
    def add(self,coeffs):
        # calculate pressure force
        gradPw = (self._P[1:-1] - self._P[0:-2])/self._grid.dx_WP
        gradPe = (self._P[2:] - self._P[1:-1])/self._grid.dx_PE
        force = 0.5*(gradPw + gradPe)*self._grid.vol

        # calculate linearization coefficients
        coeffW = -0.5*self._grid.vol/self._grid.dx_WP
        coeffE = 0.5*self._grid.vol/self._grid.dx_PE
        coeffP = -coeffW - coeffE


        # modify the linearization coefficient on the boundaries
        coeffP[0] += coeffW[0]*self._west_bc.coeff()
        coeffP[-1] += coeffE[-1]*self._east_bc.coeff()


        # zero the boundary coefficient that are not used
        coeffW[0] = 0.0
        coeffE[-1] = 0.0

        # add to coefficient arrays
        coeffs.accumulate_aP(coeffP)
        coeffs.accumulate_aW(coeffW)
        coeffs.accumulate_aE(coeffE)
        coeffs.accumulate_rP(force)

        # return the modified coefficient array
        return coeffs

# class define a mass conservation equation
class MassConservationEquation:
    # constructor
    def __init__(self,grid,U,P,dhat,Uhe,rho,P_west_bc,P_east_bc, U_west_bc,U_east_bc):
        self._grid = grid
        self._U = U
        self._P = P
        self._dhat = dhat
        self._Uhe = Uhe
        self._rho = rho
        self._P_west_bc = P_west_bc
        self._P_east_bc = P_east_bc
        self._U_west_bc = U_west_bc
        self._U_east_bc = U_east_bc

    # add diffusion term to coefficient arrays
    def add(self,PP_coeffs,PU_coeffs):
        # calculate mass imbalance, using advecting velocities
        imbalance = self._rho*self._grid.Ae*self._Uhe[1:] - self._rho*self._grid.Aw*self._Uhe[:-1]

        # calculate linearization coefficient on pressure
        PP_coeffW = np.concatenate((np.array([0]),-self._rho*self._grid.Aw[1:]*self._dhat[1:-1]/
                                    self._grid.dx_WP[1:]))
        
        PP_coeffE = np.concatenate((-self._rho*self._grid.Ae[:-1]*self._dhat[1:-1]/
                                    self._grid.dx_PE[:-1],np.array([0])))
        
        PP_coeffP = -PP_coeffW - PP_coeffE

        # calculate the linearization coefficients on velocity
        PU_coeffW = np.concatenate((np.array([-self._rho*self._grid.Aw[0]]),
                                     -0.5*self._rho*self._grid.Aw[1:]))

        PU_coeffE = np.concatenate((0.5*self._rho*self._grid.Ae[:-1],
                                     np.array([self._rho*self._grid.Ae[-1]])))

        PU_coeffP = np.concatenate((np.array([0]),PU_coeffW[1:]))+ np.concatenate((PU_coeffE[:-1],np.array([0])))

        # modify linearization coefficient on the boundaries
        # (only velocity, because pressure is already zero)
        PU_coeffP[0] += PU_coeffW[0]*self._U_west_bc.coeff()
        PU_coeffP[-1] += PU_coeffE[-1]*self._U_east_bc.coeff()

        # zero the boundary coefficients that are not used
        PU_coeffW[0] = 0.0
        PU_coeffE[-1] = 0.0

        # add to coefficient arrays
        PP_coeffs.accumulate_aP(PP_coeffP)
        PP_coeffs.accumulate_aW(PP_coeffW)
        PP_coeffs.accumulate_aE(PP_coeffE)
        PP_coeffs.accumulate_rP(imbalance)

        PU_coeffs.accumulate_aP(PU_coeffP)
        PU_coeffs.accumulate_aW(PU_coeffW)
        PU_coeffs.accumulate_aE(PP_coeffE)

        # return the modified coefficient arrays
        return PP_coeffs, PU_coeffs

# class define an extrapolated boundary conditions
# class ExtrapolateBC:
#     # constructor
#     def __init__(self,phi,grid,loc):
#         self._phi = phi
#         self._grid = grid
#         self._loc = loc

#     # return boundary condition value
#     def value(self):
#         if self._loc is BoundaryLocation.WEST:
            
        




## MAIN CODE ##
lx = 4.0
ly = 0.02
lz = 0.02
ncv = 10
grid = Grid(lx,ly,lz,ncv)

# set timestep
nTime = 1
dt = 1e9
time = 0

# iteration info
maxIter = 100
converged = 1e-6

# thermo properties
rho = 1000
mu = 1e-3

# coefficients
PU_coeffs = ScalarCoeffs(grid.ncv)
PP_coeffs = ScalarCoeffs(grid.ncv)
UP_coeffs = ScalarCoeffs(grid.ncv)
UU_coeffs = ScalarCoeffs(grid.ncv)

# initial condition
U0 = 10
P0 = 0

# initialize field variable arrays
U = U0*np.ones(grid._ncv+2)
P = P0*np.ones(grid._ncv+2)

# initialize advecting and damping coefficient arrays
dhat = np.zeros(grid.ncv+1)
Uhe = np.zeros(grid.ncv+1)

# define boundary conditions for velocity
U_west_bc = DirichletBC(U,grid,U0,BoundaryLocation.WEST)
U_east_bc = NeumannBC(U,grid,0,BoundaryLocation.EAST)


# define boundary conditions for pressure
P_west_bc = NeumannBC(P,grid,0,BoundaryLocation.WEST)
P_east_bc = DirichletBC(P,grid,0,BoundaryLocation.EAST)


# apply boundary conditions
U_west_bc.apply()
U_east_bc.apply()
P_west_bc.apply()
P_east_bc.apply()

# define the transient mdel
Uold = np.copy(U)
transient = FirstOrderTransientModel(grid,U,Uold,rho,1,dt)

# define the diffusion model
diffusion = DiffusionModel(grid,U,mu,U_west_bc,U_east_bc)

#define advection model
advection = UpwindAdvectionModel(grid,U,Uhe,rho,1,U_west_bc,U_east_bc)

# define the pressure force model
pressure = PressureForceModel(grid,P, P_west_bc, P_east_bc)

# define the advecting velocity model
advecting = AdvectingVelocityModel(grid,dhat,Uhe,P,U,UU_coeffs)

# define conservation of mass equation
mass = MassConservationEquation(grid,U,P,dhat,Uhe,rho,P_west_bc,
                                P_east_bc,U_west_bc,U_east_bc)

# Loop through all timesteps
for tStep in range(nTime):
    time += dt
    print("Timestep = {}; Time = {}".format(tStep,time))

    Uold[:] = U[:]

    for i in range(maxIter):
        # zero out coefficients
        PP_coeffs.zero()
        PU_coeffs.zero()
        UU_coeffs.zero()
        UP_coeffs.zero()

        # assemble the momentum equation
        UU_coeffs = diffusion.add(UU_coeffs)
        UU_coeffs = advection.add(UU_coeffs)
        UU_coeffs = transient.add(UU_coeffs)
        UU_coeffs = pressure.add(UU_coeffs)

        # assemble the mass equation
        advecting.update()
        PP_coeffs,PU_coeffs = mass.add(PP_coeffs,PU_coeffs)

        # compute the residual and check for convergence
        PmaxResid = norm(PU_coeffs.rP + PP_coeffs.rP,np.inf)
        PavgResid = np.mean(np.absolute(PU_coeffs.rP + PP_coeffs.rP))
        
        UmaxResid = norm(UU_coeffs.rP + UU_coeffs.rP,np.inf)
        UavgResid = np.mean(np.absolute(UU_coeffs.rP + UU_coeffs.rP))

        print("Iteration = {}".format(i))
        print("Mass:  Max.Resid = {}; Avg.Resid = {}".format(PmaxResid,PavgResid))
        print("Momentum:  Max.Resid = {}; Avg.Resid = {}".format(UmaxResid,UavgResid))

        if PmaxResid < converged and UmaxResid < converged:
            break

        # solve the sparse matrix system
        dP, dU = solve(PP_coeffs,PU_coeffs, UP_coeffs, UU_coeffs)

        # update the solution
        P[1:-1] += dP
        U[1:-1] += dU

        # update boundary condition
        U_west_bc.apply()
        U_east_bc.apply()
        P_west_bc.apply()
        P_east_bc.apply()


        # update the advectig velocities
        advecting.update()
        

print(np.shape(U))
print(np.shape(P))
print(np.shape(grid.xP))

plt.plot(grid.xP, U)
plt.show()
        
        
    



        
        
