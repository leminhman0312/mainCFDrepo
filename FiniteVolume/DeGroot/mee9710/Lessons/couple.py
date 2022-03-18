from Lesson5.Grid import Grid
from Lesson5.ScalarCoeffs import ScalarCoeffs
from Lesson5.BoundaryConditions import BoundaryLocation, DirichletBc, NeumannBc
from Lesson5.Models import DiffusionModel
from Lesson5.LinearSolver import solve
from Lesson5.Models import FirstOrderTransientModel, DiffusionModel, UpwindAdvectionModel, PressureForceModel
from Lesson5.Models import AdvectingVelocityModel, MassConservationEquation

import numpy as np
from numpy.linalg import norm


## MAIN CODE ##
lx = 4.0
ly = 0.02
lz = 0.02
ncv = 10
grid = Grid(lx,ly,lz,ncv)

# set timestep
nTime = 1
dt = 1e-2
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
U_west_bc = DirichletBc(U,grid,U0,BoundaryLocation.WEST)
U_east_bc = NeumannBc(U,grid,0,BoundaryLocation.EAST)


# define boundary conditions for pressure
P_west_bc = NeumannBc(P,grid,0,BoundaryLocation.WEST)
P_east_bc = DirichletBc(P,grid,0,BoundaryLocation.EAST)


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
        UP_coeffs = pressure.add(UP_coeffs)

        # assemble the mass equation
        advecting.update()
        PP_coeffs,PU_coeffs = mass.add(PP_coeffs,PU_coeffs)

        # compute the residual and check for convergence
        PmaxResid = norm(PU_coeffs.rP + PP_coeffs.rP,np.inf)
        PavgResid = np.mean(np.absolute(PU_coeffs.rP + PP_coeffs.rP))
        
        UmaxResid = norm(UU_coeffs.rP + UP_coeffs.rP,np.inf)
        UavgResid = np.mean(np.absolute(UU_coeffs.rP + UP_coeffs.rP))

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


import matplotlib.pyplot as plt

plt.plot(grid.xP,U)
plt.show()
