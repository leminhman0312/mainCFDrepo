import numpy as np

class DiffusionModel:
    """Class defining a diffusion model"""

    def __init__(self, grid, phi, gamma, west_bc, east_bc):
        """Constructor"""
        self._grid = grid
        self._phi = phi
        self._gamma = gamma
        self._west_bc = west_bc
        self._east_bc = east_bc

    def add(self, coeffs):
        """Function to add diffusion terms to coefficient arrays"""

        # Calculate the west and east face diffusion flux terms for each face
        flux_w = - self._gamma*self._grid.Aw*(self._phi[1:-1]-self._phi[0:-2])/self._grid.dx_WP
        flux_e = - self._gamma*self._grid.Ae*(self._phi[2:]-self._phi[1:-1])/self._grid.dx_PE

        # Calculate the linearization coefficients
        coeffW = - self._gamma*self._grid.Aw/self._grid.dx_WP
        coeffE = - self._gamma*self._grid.Ae/self._grid.dx_PE
        coeffP = - coeffW - coeffE

        # Modify the linearization coefficients on the boundaries
        coeffP[0] += coeffW[0]*self._west_bc.coeff()
        coeffP[-1] += coeffE[-1]*self._east_bc.coeff()

        # Zero the boundary coefficients that are not used
        coeffW[0] = 0.0
        coeffE[-1] = 0.0

        # Calculate the net flux from each cell
        flux = flux_e - flux_w

        # Add to coefficient arrays
        coeffs.accumulate_aP(coeffP)
        coeffs.accumulate_aW(coeffW)
        coeffs.accumulate_aE(coeffE)
        coeffs.accumulate_rP(flux)

        # Return the modified coefficient array
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
