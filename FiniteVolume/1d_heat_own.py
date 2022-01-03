# Solving the 1D Diffusion Equation
import numpy as np
from enum import Enum
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from numpy.linalg import norm
import sys
import matplotlib.pyplot as plt


def accumulate_aP(aP):
    # accumulate values onto aP
    aP += aP

def accumulate_aW(aW):
    # accumulate values onto aW
    aW += aW

def accumulate_aE(aE):
    # accumulate values onto aE
    aE += aE

def accumulate_rP(rP):
    # accumulate values onto rP
    rP += rP

def DirichletBC_apply(phi,value,loc):
    if loc == 1:                # west
        phi[0] = value
    elif loc == 2:
        phi[-1] = value
    else:
        raise ValueError("Unknown boundary location")

    return 1

        

def DiffusionAdd(phi,gamma, dx_WP, dx_PE, bc_left, bc_right, Aw, Ae, rP):
    # calculate west/east diffusion flux for each face
    flux_w = -gamma*Aw*(phi[1:-1]-phi[0:-2])/dx_WP 
    flux_e = -gamma*Ae*(phi[2:]-phi[1:-1])/dx_PE
    
    # calculate the linearized coefficient [aW, aP, aE]
    coeffW = -gamma*Aw/dx_WP
    print('coeffW = ', coeffW,'gamma = ', gamma, ' dx_WP = ', dx_WP)
    coeffE = -gamma*Ae/dx_PE
    coeffP = -coeffW - coeffE
    
    # modified the linearized coefficient on the boundaries
    coeffP[0] += coeffW[0]*bc_left
    coeffE[-1] += coeffE[-1]*bc_right

    # zero the coefficients that are not used (on the boundary)
    coeffW[0] = 0.0
    coeffE[-1] = 0.0

    # calculate the net flux from each cell (out -in)
    flux = flux_e - flux_w

    # add to coefficient arrays
    accumulate_aP(coeffP)
    accumulate_aW(coeffW)
    accumulate_aE(coeffE)
    accumulate_rP(flux)
   


# function to return a sparse matrix representation of a set of scalar coefficients
def get_sparse_matrix(ncv, aE, aP, aW):
                      
    ncv = ncv
    data = np.zeros(3*ncv-2)
    rows = np.zeros(3*ncv-2, dtype = int)
    cols = np.zeros(3*ncv-2, dtype = int)
    data[0] = aP[0]
    rows[0] = 0
    cols[0] = 0

    if ncv > 1:
        data[1] = aE[0]
        rows[1] = 0
        cols[1] = 1

    for i in range(ncv-2):
        data[3*i+2] = aW[i+1]
        data[3*i+3] = aP[i+1]
        data[3*i+4] = aE[i+1]

        rows[3*i+2:3*i+5] = i+1

        cols[3*i+2] = i
        cols[3*i+3] = i+1
        cols[3*i+4] = i+2

    if ncv > 1:
        data[3*ncv-4] = aW[-1]
        data[3*ncv-3] = aP[-1]
        
        rows[3*ncv-4:3*ncv-2] = ncv-1

        cols[3*ncv-4] = ncv-2
        cols[3*ncv-3] = ncv-1

    return csr_matrix((data, (rows,cols)))

def solve(A, rP):
    return spsolve(A, -rP)
    




## MAIN CODE ##
# Solving the 1D steady conduction with Dirichet BC
# Initially, east = 300K, west = 300K
# define the grid
lx = 1.0
ly = 0.1
lz = 0.1
ncv = 10
# calculate the control volume length
dx = lx/float(ncv)

# calculate the face locations
# generate array of control volume element with length dx each
xf = np.array([i*dx for i in range(ncv+1)])

# calculate the cell centroid locations
# Note: figure out why add a xf[-1] last item 
xP = np.array([xf[0] +
                        [0.5*(xf[i]+xf[i+1]) for i in range(ncv)] +
                        [xf[-1]]])
# calculate face areas
Af = ly*lz*np.ones(ncv+1)

# calculate the outer surface area of each cell
Ao = (2.0*dx*ly + 2.0*dx*lz)*np.ones(ncv)

# calculate cell volumes
vol = dx*ly*lz*np.ones(ncv)

dx_WP = dx*np.ones(ncv) #xP[1:-1]-xP[0:-2]


dx_PE = dx*np.ones(ncv) # xP[2:]-xP[1:-1]



Aw = Af[0:-1]

Ae = Af[1:]

#set the max iterations and convergence criterion
maxIter = 100
converged=1e-6

# thermal properties
k = 0.1


#create ap, aw, ae, rp array
ncv = ncv
aP = np.zeros(ncv)
aW = np.zeros(ncv)
aE = np.zeros(ncv)
rP = np.zeros(ncv)

# zero them out
aP.fill(0.0)
aW.fill(0.0)
aE.fill(0.0)
rP.fill(0.0)

# initial condition
T0 = 300

# initialize field variable arrays
T = T0*np.ones(ncv+2)

# set the boundary conditions
left = DirichletBC_apply(T, 400, 1)    # west
right = DirichletBC_apply(T, 300, 2)    # east


# array to store solutions at each iterations
T_solns = [np.copy(T)]


# # iterate until converged
for i in range(maxIter):
    # print ('i = ', i)
    # zero out coefficients
    aP.fill(0.0)
    aW.fill(0.0)
    aE.fill(0.0)
    rP.fill(0.0)

    # fill out the coefficient
    DiffusionAdd(T,k, dx_WP, dx_PE, left, right, Aw, Ae , rP)

    # compute residual and check for convergence
    maxResid = norm(rP, np.inf)
    avgResid = np.mean(np.absolute(rP))
    # print("Iteration = {} ; Resid = {} ; Avg. Resid = {}".format(i,maxResid, avgResid) )
    if maxResid > converged:
        break


    # get the matrix
    A = get_sparse_matrix(ncv, aE, aP, aW)
    # solve the sparse matrix system
    dT = solve(A, rP)


    # update the solution
    T[1:-1] += dT
    left = DirichletBC_apply(T, 400, 1)    # west
    right = DirichletBC_apply(T, 300, 2)    # east

    # store the solution
    T_solns.append(np.copy(T))




# plotting
# i = 0
# print(T_solns)


# plt.plot(xP[0],T_solns[0:11])
# plt.xlabel("X")
# plt.ylabel("T")
# plt.legend()
# plt.show()
    
