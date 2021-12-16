!! Solution of 2-D Euler- and Navier-Stokes Equations
!! on Unstructured Triangular Grids.
!
!  Features:
!  ~~~~~~~~~
!  # unstructured finite-volume scheme of median-dual type
!  # triangular elements only
!  # ideal gas model (other models possible)
!  # laminar flow (viscosity computed by Sutherland`s law)
!  # Roe's flux-difference splitting scheme, Venkatakrishnan's limiter
!  # explicit multistage time-stepping scheme (Runge-Kutta)
!  # preconditioning for low Mach numbers
!  # global or local time steps
!  # central implicit residual smoothing
!  # characteristic boundary conditions for external and internal flows
!  # special initial solution for compressor and turbine blades
!
!  (c) J. Blazek, CFD Consulting & Analysis, www.cfd-ca.de
!  Created February 25, 2014
!  Version 1.0 from September 19, 2014
!
! *****************************************************************************
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!
! *****************************************************************************
module ModDataTypes
  implicit none
  integer, parameter :: chrlen = 256         !< length of strings
  integer, parameter :: rtype  = kind(1.D0)  !< reals with double precision
end module ModDataTypes
module ModControl
  use ModDataTypes
  implicit none
  character(1) :: lrest !< use of previous solution for restart ("Y"=yes, "N"=no)
  integer :: maxiter, & !< max. number of iterations
             outstep, & !< number of iterations between solution dumps
             iter       !< actual iteration number
  real(rtype) :: convtol !< convergence criterion (2-norm of density change for
                         !! which the iteration process is stopped)
end module ModControl
module ModFiles
  use ModDataTypes
  implicit none
  character(chrlen) :: fnGrid, & !< grid and topology data
                       fnFlow, & !< flow field (+ 5 digit iteration number + .v2d)
                       fnSurf, & !< quantities along wall surface(s) (+ 5 digit iteration number + .v2d)
                       fnConv, & !< convergence history (+ .v2d)
                       fnRsti, & !< restart solution - input
                       fnRsto    !< restart solution - output
  integer, parameter :: ifInp  = 10, & !< user input file (name stored in main.f90)
                        ifGrid = 20, &
                        ifFlow = 30, &
                        ifSurf = 40, &
                        ifConv = 50, &
                        ifRsti = 60, &
                        ifRsto = 70
end module ModFiles
module ModGeometry
  use ModDataTypes
  implicit none
  integer :: nnodes !< number of grid nodes (including dummy nodes)\n
    !! @details
    !! Dummy nodes are defined for inlet, outlet and far-field boundaries only. The
    !! variables cv(), dv(), x(), y(), edge() and sij() are divided into two parts.
    !! The first part is related to the physical grid (dimensions nndint, nedint). The
    !! second part refers to the dummy nodes. This approach makes it easier to loop
    !! either only over the physical grid nodes (edges), or over all nodes (edges)
    !! using the same vector.
  integer :: nndint, &  !< number of physical grid nodes (nnodes - dummy nodes)
             ntria, &   !< number of triangles
             nsegs, &   !< number of boundary segments (composed of boundary faces)
             nbfaces, & !< number of all boundary faces
             nbnodes    !< total number of boundary nodes
  character(chrlen), allocatable :: bname(:) !< names of boundary segments (used for plots)
  integer, allocatable :: btype(:) !< types of boundary conditions:\n
    !! 100-199 = inflow\n
    !! 200-299 = outflow\n
    !! 300-399 = viscous wall\n
    !! 400-499 = inviscid wall\n
    !! 500-599 = symmetry line\n
    !! 600-699 = far-field\n
    !! 700-799 = periodic boundary
  integer, allocatable :: bface(:,:) !< indexes of two nodes defining a face (NOT used for periodic boundaries!)
  integer, allocatable :: bnode(:,:) !< data related to boundary nodes\n
    !! @details
    !! Meaning of the entries:
    !! @li (1,*) = index of the node itself
    !! @li (2,*) = index of the related dummy node (inlet, outlet, far-field)
    !!             or of the other periodic node; otherwise = -777
    !! @li (3,*) = index of edge to the dummy node (edge(*,ie), ie > nedint)
  integer, allocatable :: ibound(:,:) !< pointer from boundary segment to boundary faces and nodes\n
    !! @details
    !! Meaning of the entries:
    !! @li (1,*) = last index in bface()
    !! @li (2,*) = last index in bnode()
  real(rtype) :: xref, & !< x-coordinate of the reference point
                 yref, & !< y-coordinate of the reference point
                 cref    !< reference length or airfoil chord
  integer, allocatable :: tria(:,:) !< node indexes of triangle elements
  real(rtype), allocatable :: x(:), &  !< x-coordinates of grid points
                              y(:)     !< y-coordinates of grid points
  real(rtype), allocatable :: sij(:,:) !< x,y-components of the face vector (n*dS)\n
    !! @details
    !! For ie > nedint, sij(*,ie) represents the average face vector at a
    !! boundary node (face between boundary and dummy node);
    !! sij() always points from node i to node j (see Fig. 5.9).
  real(rtype), allocatable :: vol(:)   !< median-dual control volume (shaded area in Fig. 5.8)
  real(rtype), allocatable :: sbf(:,:) !< normal vector of boundary face (outward pointing,
    !! size equal to the edge length); defined for all boundaries except for periodic ones
  real(rtype), allocatable :: sproj(:,:) !< projections of control volumes on the x- and y-axis
    !! (used to compute the time step - see Eq. (6.22))
end module ModGeometry
module ModNumerics
  use ModDataTypes
  implicit none
  character(1) :: ktimst,  & !< switch between local (="L") and global (="G") time-stepping
                  lvort,   & !< far-field vortex correction ("Y"=yes, "N"=no)
                  kprecond   !< low Mach-number preconditioning ("Y"=yes, "N"=no)
  integer :: nedges, & !< total number of edges (including edges between boundary and dummy nodes)
             nedint, & !< number of edges excluding those to dummy nodes
             iorder, & !< order of Roe's upwind scheme (1 or 2)
             nitirs, & !< number of Jacobi iterations (implicit residual smoothing)
             nrk,    & !< number of stages (Runge-Kutta scheme); max. = 5
             ldiss(5)  !< dissipation evaluation per stage (0=no, 1=yes)
  integer, allocatable :: edge(:,:) !< edge list (node i, node j)\n
    !! @details
    !! For ie > nedint, edge(*,ie) represents the edge from a boundary node
    !! (i) to a dummy node (used at inlet, outlet and far-field boundaries).
  real(rtype) :: cfl,      & !< CFL-number
                 epsirs,   & !< coefficient of implicit residual smoothing
                 limfac,   & !< limiter coefficient (Roe's upwind scheme)
                 epsentr,  & !< entropy correction coefficient (Roe's upwind scheme)
                 precoeff, & !< preconditioning parameter K (low Mach numbers)
                 ark(5),   & !< stage coefficients
                 betrk(5)    !< dissipation-blending coefficients
  real(rtype) :: volref    !< reference volume\n
                           !! @details
                           !! Parameter is required for the computation of limiter
                           !! functions (higher-order Roe scheme).
  real(rtype) :: limref(4) !< reference values of density, u, v and pressure\n
                           !! @details
                           !! Parameter is required for the computation of limiter
                           !! functions (higher-order Roe scheme).
  real(rtype), allocatable :: cvold(:,:), & !< conservative variables from previous time step
                              diss(:,:),  & !< artificial dissipation
                              rhs(:,:),   & !< residual (right-hand side)
                              lim(:,:),   & !< values of the limiter function (density, u, v, pressure)
                              tstep(:)      !< time steps (without the CFL-number)
  real(rtype), allocatable :: gradx(:,:) !< gradients of density, velocity components,
                              !! pressure and temperature with respect to the x-coordinate
  real(rtype), allocatable :: grady(:,:) !< gradients of density, velocity components,
                              !! pressure and temperature with respect to the y-coordinate
  real(rtype) :: pi, & !< 3.14...
                 rad   !< 180./pi
end module ModNumerics
module ModPhysics
  use ModDataTypes
  implicit none
  character(1) :: kequs, & !< equations solved ("E"=Euler, "N"=Navier-Stokes)
                  kflow    !< type of flow ("E"=external, "I"=internal)
! reference values
  real(rtype) :: gamma,  & !< ratio of specific heat coefficients
                 cpgas,  & !< specific heat coefficient at constant pressure
                 prlam,  & !< laminar Prandtl number
                 renum,  & !< Reynolds number
                 refvel, & !< reference velocity (internal flow only; for external flow computed from the far-field boundary)
                 refrho, & !< reference density (internal flow only; for external flow computed from the far-field boundary)
                 refvisc   !< reference dynamic viscosity coefficient (computed from renum, refvel, cref and refrho)
! boundary conditions - external flow
  real(rtype) :: machinf, & !< Mach-number at infinity
                 alpha,   & !< angle of attack
                 pinf,    & !< static pressure at infinity
                 tinf,    & !< static temperature at infinity
                 rhoinf,  & !< density at infinity
                 uinf,    & !< u-component of velocity vector at infinity
                 vinf,    & !< v-component of velocity vector at infinity
                 qinf       !< total velocity (= SQRT(uinf**2+vinf**2))
! boundary conditions - internal flow
  real(rtype) :: ptinl,   & !< total pressure at inlet
                 ttinl,   & !< total temperature at inlet
                 betainl, & !< low angle at inlet (with x-axis, positive in the clock-wise direction)
                 betaout, & !< approximate outlet angle (utilized for the initial guess only)
                 p12rat,  & !< ratio of inlet to outlet static pressure (initial guess only)
                 pout       !< static pressure at outlet
! flow variables
  integer :: nconv, & !< number of conservative variables (cv)
             ndepv    !< number of dependent variables (dv)
  real(rtype), allocatable :: cv(:,:) !< conservative variables\n
              !! @details
              !! cv(1,i) = density\n
              !! cv(2,i) = density * u\n
              !! cv(3,i) = density * v\n
              !! cv(4,i) = density * E
  real(rtype), allocatable :: dv(:,:) !< dependent variables\n
              !! @details
              !! dv(1,i) = static pressure\n
              !! dv(2,i) = static temperature\n
              !! dv(3,i) = speed of sound\n
              !! dv(4,i) = ratio of specific heats\n
              !! dv(5,i) = specific heat coefficient at constant pressure\n
              !! dv(6,i) = laminar viscosity coefficient (if viscous flow)\n
              !! dv(7,i) = laminar heat conductivity coefficient (if viscous flow)
end module ModPhysics
module ModPlotQuant
  use ModDataTypes
  implicit none
  integer, parameter :: mxquant =13, & !< total number of plot variables
                        mxqfield=11    !< no. of plot variables in the field (cf and Cp only at the boundaries)
  character(chrlen) :: title           !< title of the simulation case
  character(chrlen) :: cquant(mxquant) !< names of plot variables
  character(1)      :: lquant(mxquant) !< on/off switches of the plot variables
  real(rtype) :: drho,  & !< change of the density residual (convergence criterion)
                 drho1, & !< initial change of the density residual (used for normalization)
                 cl,    & !< lift coefficient (pressure forces only; external flow)
                 cd,    & !< drag coefficient (pressure forces only; external flow)
                 cm,    & !< pitching moment coefficient wrp. to the reference point
                 mflow, & !< average mass flow rate (internal flow)
                 mfratio  !< ratio of mass flow at outlet to mass flow at inlet
end module ModPlotQuant
module ModInterfaces
  implicit none
  interface
  subroutine AllocateMemory
  end subroutine AllocateMemory
  subroutine BcondFarfield( ibegn,iendn,rhof,uf,vf,pf)
    use ModDataTypes
    integer, intent(in) :: ibegn, iendn
    real(rtype) :: rhof(:), uf(:), vf(:), pf(:)
  end subroutine BcondFarfield
  subroutine BcondInflow( ibegn,iendn)
    integer, intent(in) :: ibegn, iendn
  end subroutine BcondInflow
  subroutine BcondOutflow( ibegn,iendn)
    integer, intent(in) :: ibegn, iendn
  end subroutine BcondOutflow
  subroutine BcondWallns( ibegn,iendn)
    integer, intent(in) :: ibegn, iendn
  end subroutine BcondWallns
  subroutine BoundaryConditions( work)
    use ModDataTypes
    real(rtype) :: work(:)
  end subroutine BoundaryConditions
!djb
! subroutine CheckMetrics( work)
!   use ModDataTypes
!   real(rtype) :: work(:)
! end subroutine CheckMetrics
!djb
  function CompTheta( gam,c,q2)
    use ModDataTypes
    real(rtype), intent(in) :: gam, c, q2
    real(rtype) :: CompTheta
  end function CompTheta
  subroutine Cons2Prim( wvec,wpvec,H,q2,theta,rhoT,hp,hT,pmat1)
    use ModDataTypes
    real(rtype), intent( in) :: H, q2, theta, rhoT, hp, hT
    real(rtype), intent( in) :: wvec(5), wpvec(5)
    real(rtype), intent(out) :: pmat1(5,5)
  end subroutine Cons2Prim
  subroutine Convergence
  end subroutine Convergence
  subroutine DependentVarsAll
  end subroutine DependentVarsAll
  subroutine DependentVarsOne( i)
    integer, intent(in) :: i
  end subroutine DependentVarsOne
  subroutine DissipRoe1( beta)
    use ModDataTypes
    real(rtype), intent(in) :: beta
  end subroutine DissipRoe1
  subroutine DissipRoe1Prec( beta)
    use ModDataTypes
    real(rtype), intent(in) :: beta
  end subroutine DissipRoe1Prec
  subroutine DissipRoe2( beta)
    use ModDataTypes
    real(rtype), intent(in) :: beta
  end subroutine DissipRoe2
  subroutine DissipRoe2Prec( beta)
    use ModDataTypes
    real(rtype), intent(in) :: beta
  end subroutine DissipRoe2Prec
  subroutine DummyNodes
  end subroutine DummyNodes
  subroutine EdgesFinalize( niedge,iedge)
    integer :: niedge(:), iedge(:,:)
  end subroutine EdgesFinalize
  subroutine EdgesInitialize( niedge,iedge)
    integer, intent(out) :: niedge(:), iedge(:,:)
  end subroutine EdgesInitialize
  subroutine ErrorMessage( message)
    character(*), intent(in) :: message
  end subroutine ErrorMessage
  subroutine FaceVectorsSymm( marker)
    integer :: marker(:)
  end subroutine FaceVectorsSymm
  subroutine FluxViscous( beta)
    use ModDataTypes
    real(rtype), intent(in) :: beta
  end subroutine FluxViscous
  subroutine FluxWalls
  end subroutine FluxWalls
  subroutine Forces
  end subroutine Forces
  subroutine FluxRoe1
  end subroutine FluxRoe1
  subroutine FluxRoe2
  end subroutine FluxRoe2
  subroutine Gradients
  end subroutine Gradients
  subroutine GradientsVisc
  end subroutine GradientsVisc
  subroutine InitConstants
  end subroutine InitConstants
  subroutine InitMetrics( niedge,iedge)
    integer, intent(in) :: niedge(:), iedge(:,:)
  end subroutine InitMetrics
  subroutine InitMetricsBound( marker,btria)
    integer :: marker(:), btria(:,:)
  end subroutine InitMetricsBound
  subroutine InitSolution
  end subroutine InitSolution
  subroutine Irsmoo( ncontr,rhsold,rhsit)
    use ModDataTypes
    integer     :: ncontr(:)
    real(rtype) :: rhsold(:,:), rhsit(:,:)
  end subroutine Irsmoo
  subroutine LeftEigenvec( wvec,wpvec,nvec,V,theta,rhop,rhoT,hp,hT,evl)
    use ModDataTypes
    real(rtype), intent( in) :: V, theta, rhop, rhoT, hp, hT
    real(rtype), intent( in) :: wvec(5), wpvec(5), nvec(3)
    real(rtype), intent(out) :: evl(5,5)
  end subroutine LeftEigenvec
  subroutine Limiter( umin,umax)
    use ModDataTypes
    real(rtype), intent(in) :: umin(:,:), umax(:,:)
  end subroutine Limiter
  subroutine LimiterInit( umin,umax)
    use ModDataTypes
    real(rtype), intent(out) :: umin(:,:), umax(:,:)
  end subroutine LimiterInit
  subroutine LimiterRefvals
  end subroutine LimiterRefvals
  subroutine Massflow
  end subroutine Massflow
  subroutine MatprodTp1_P1( wvec,wpvec,nvec,V,H,theta,rhop,rhoT,hp,hT,q2,mat)
    use ModDataTypes
    real(rtype), intent( in) :: V, H, theta, rhop, rhoT, hp, hT, q2
    real(rtype), intent( in) :: wvec(5), wpvec(5), nvec(3)
    real(rtype), intent(out) :: mat(5,5)
  end subroutine MatprodTp1_P1
  subroutine MatrixTimesInverse( wpvec,q2,amat,bmat,cmat)
    use ModDataTypes
    real(rtype), intent( in) :: q2
    real(rtype), intent( in) :: wpvec(5), amat(5,5), bmat(5,5)
    real(rtype), intent(out) :: cmat(5,5)
  end subroutine MatrixTimesInverse
  subroutine MatVecProd5( a,v,c)
    use ModDataTypes
    real(rtype), intent( in) :: a(5,5), v(5)
    real(rtype), intent(out) :: c(5)
  end subroutine MatVecProd5
  subroutine Periodic( var)
    use ModDataTypes
    real(rtype) :: var(:,:)
  end subroutine Periodic
  subroutine PlotFlow
  end subroutine PlotFlow
  subroutine PlotSurfaces
  end subroutine PlotSurfaces
  subroutine Prim2Cons( wvec,wpvec,H,theta,rhoT,hp,hT,pmat)
    use ModDataTypes
    real(rtype), intent( in) :: H, theta, rhoT, hp, hT
    real(rtype), intent( in) :: wvec(5), wpvec(5)
    real(rtype), intent(out) :: pmat(5,5)
  end subroutine Prim2Cons
  subroutine PrintParams
  end subroutine PrintParams
  function ReadChar( iunit)
    integer, intent(in) :: iunit
    character(1) :: ReadChar
  end function ReadChar
  subroutine ReadGrid
  end subroutine ReadGrid
  subroutine ReadParams( fname)
    character(*), intent(in) :: fname
  end subroutine ReadParams
  subroutine ReadSolution
  end subroutine ReadSolution
  subroutine RightEigenvec( wvec,wpvec,nvec,V,H,theta,rhop,rhoT,hp,hT,evr)
    use ModDataTypes
    real(rtype), intent( in) :: V, H, theta, rhop, rhoT, hp, hT
    real(rtype), intent( in) :: wvec(5), wpvec(5), nvec(3)
    real(rtype), intent(out) :: evr(5,5)
  end subroutine RightEigenvec
  subroutine Solver( iwork,work)
    use ModDataTypes
    integer     :: iwork(:)
    real(rtype) :: work(:)
  end subroutine Solver
  subroutine TimeStep
  end subroutine TimeStep
  subroutine VolumeProjections
  end subroutine VolumeProjections
  subroutine WriteSolution
  end subroutine WriteSolution
  subroutine ZeroResiduals
  end subroutine ZeroResiduals
  end interface
end module ModInterfaces
!> Main program of the flow solver.
!!
program Unstruct2D
  use ModDataTypes
  use ModControl
  use ModFiles
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces
  implicit none
! local variables
  character(chrlen) :: fname               ! filename (user input, convergence)
  integer           :: errFlag             ! error flag
  integer, allocatable     :: niedge(:), & ! temporary edge structures (DO NOT overwrite
                              iedge(:,:)   ! between EdgesInitialize and InitMetrics)
  integer, allocatable     :: iwork(:)     ! integer work space (used by Irsmoo)
  real(rtype), allocatable :: work(:)      ! real work space (used inside Solver)
! *****************************************************************************
  write(*,*)'*************************************************'
  write(*,*)'*                                               *'
  write(*,*)'*   2-D FLOW ON UNSTRUCTURED TRIANGULAR GRIDS   *'
  write(*,*)'*                                               *'
  write(*,*)'*  (c) Jiri Blazek, CFD Consulting & Analysis   *'
  write(*,*)'*                 www.cfd-ca.de                 *'
  write(*,*)'*                                               *'
  write(*,*)'*          Version 1.0 from 09/19/2014          *'
  write(*,*)'*                                               *'
  write(*,*)'*************************************************'
! read name of input file from command line
  call Getarg(1,fname)
  if(Len_trim(fname) == 0) call Usage
! set names of quantities which can be written out
  cquant( 1) = "density [kg/m^3]" ! density
  cquant( 2) = "u [m/s]"          ! u-velocity
  cquant( 3) = "v [m/s]"          ! v-velocity
  cquant( 4) = "p [Pa]"           ! static pressure
  cquant( 5) = "p-tot [Pa]"       ! total pressure
  cquant( 6) = "T [K]"            ! static temperature
  cquant( 7) = "T-tot [K]"        ! total temperature
  cquant( 8) = "M"                ! local Mach-number
  cquant( 9) = "M-isen"           ! isentropic Mach-number
  cquant(10) = "pt-loss"          ! total pressure loss
  cquant(11) = "visc-lam"         ! laminar viscosity
  cquant(12) = "cf"               ! friction coefficient (boundaries only)
  cquant(13) = "-Cp"              ! pressure coefficient (boundaries only)
! read input parameters
  call ReadParams(fname)
! initialize some constants
  call InitConstants
! print input parameters for checking
  call PrintParams
! set no. of equations (rho, rho*u, rho*v, rho*E, ...);
! set no. of dependent variables (p, T, c, gamma, cpgas, ...)
  nconv = 4
  ndepv = 5
  if(kequs == "N")then
    ndepv = ndepv + 2   ! laminar viscosity, heat conduction coeff.
  endif
! read grid dimensions, coordinates, and triangle nodes
  write(*,*)'Reading grid data ...'
  call ReadGrid
! generate edge list
  write(*,*)'Generating edge list ...'
  allocate( niedge(nndint),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate temporary edge pointer")
  allocate( iedge(3,2*ntria),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate temporary edge list")
  call EdgesInitialize(niedge,iedge)
  call EdgesFinalize(niedge,iedge)
! print some statistics
  write(*,1000) nndint,nnodes-nndint,ntria,nedint,nedges,nbfaces,nbnodes
! allocate work space (real)
  allocate( work(2*(nconv*nnodes)),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for work space")
! compute face vectors & cell volumes and check them;
! set coordinates of dummy nodes = those of boundary nodes;
! correct face vectors at symmetry boundaries;
! compute projections of control volumes
  write(*,*)'Computing metrics ...'
  call InitMetrics(niedge,iedge)
  call InitMetricsBound(niedge,iedge)
  call CheckMetrics(work)
  call FaceVectorsSymm(niedge)
  deallocate( iedge )
  deallocate( niedge)
  call VolumeProjections
! allocate remaining memory
  write(*,*)'Allocating remaining memory ...'
  call AllocateMemory
  allocate( iwork(nndint),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for integer works space")
! read / initialize flow field
  if(lrest == "Y")then
    write(*,*)'Reading initial solution ...'
    call ReadSolution
  else
    write(*,*)'Guessing initial solution ...'
    call InitSolution
  endif
  call DependentVarsAll
  if(lrest /= "Y")then
    call BoundaryConditions(work)
  endif
! compute limiter reference values
  call LimiterRefvals
! open file for convergence history
  write(fname,"(A)") Trim(fnConv)//".v2d"
!!open(unit=ifConv, file=fname, status="unknown", action="write", iostat=errFlag)
!!if (errFlag /= 0) call ErrorMessage("cannot open convergence file")
!!if (kflow == "E")then
!!  write(ifConv,1010) Trim(title)
!!else
!!  write(ifConv,1020) Trim(title)
!!endif
! -----------------------------------------------------------------------------
! iterate until steady state solution or max. number of iterations is reached
! -----------------------------------------------------------------------------
  if(kflow == "E")then
    write(*,1015)
  else
    write(*,1025)
  endif
  if(lrest /= "Y") iter = 0
  do
    iter = iter + 1
    call Solver(iwork,work)
    call Convergence
    if(Mod(iter,outstep) == 0)then
      write(*,*)'Writing plot files ...'
      call PlotFlow
      call PlotSurfaces
    endif
    if(iter>=maxiter .or. drho<=convtol) exit
  enddo
! -----------------------------------------------------------------------------
! close file for convergence history
!!close(unit=ifConv)
! output the results
  if(Mod(iter,outstep) /= 0)then
    write(*,*)'Writing plot files ...'
    call PlotFlow
    call PlotSurfaces
  endif
! store solution for restart
  write(*,*)'Writing solution file ...'
  call WriteSolution
  write(*,*)'Finished.'
1000 format("No. of interior nodes: ",I8,/, &
            "No. of dummy nodes   : ",I8,/, &
            "No. of grid cells    : ",I8,/, &
            "No. of interior edges: ",I8,/, &
            "Total number of edges: ",I8,/, &
            "No. of boundary faces: ",I8,/, &
            "No. of boundary nodes: ",I8)
1010  format(A,/,"1",/,"Convergence History",/,"1 7",/, &
             "step",/,"resid",/,"resmax",/,"i-res",/,"cl",/, &
             "cd",/,"cm",/,"-1 0",/,"0 0 0",/,"Unstructured")
1015  format(75("-"),/, &
             " step",5X,"resid",7X,"resmax",5X,"i-res",6X,"cl", &
             10X,"cd",10X,"cm",/,75("-"))
1020  format(A,/,"1",/,"Convergence History",/,"1 6",/, &
             "step",/,"resid",/,"resmax",/,"i-res",/,"mass-flow",/, &
             "mass-ratio",/,"-1 0",/,"0 0 0",/,"Unstructured")
1025  format(64("-"),/, &
             " step",5X,"resid",7X,"resmax",5X,"i-res",3X,"mass flow", &
             3X,"mass ratio",/,64("-"))
! *****************************************************************************
contains
  subroutine Usage
    write(*,*)'Usage:'
    write(*,*)'Unstruct2D <input file>'
    stop
  end subroutine Usage
end program Unstruct2D
subroutine ZeroResiduals
  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none
! local variables
  integer :: i, ib, ibn, idir, ibegn, iendn
! *****************************************************************************
  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
! - symmetry boundary
    if(btype(ib)>=500 .and. btype(ib)<600)then
      if(btype(ib)-500 <  2) idir = 2  ! x=const. line -> x-component
      if(btype(ib)-500 >= 2) idir = 3  ! y=const. line -> y-component
      do ibn=ibegn,iendn
        i           = bnode(1,ibn)
        rhs(idir,i) = 0D0
      enddo
! - viscous (no-slip) wall
    else if((btype(ib)>=300 .and. btype(ib)<400) .and. kequs=="N")then
      do ibn=ibegn,iendn
        i        = bnode(1,ibn)
        rhs(2,i) = 0D0       ! velocity components = 0
        rhs(3,i) = 0D0
      enddo
    endif
    ibegn = iendn + 1
  enddo
end subroutine ZeroResiduals
subroutine FluxViscous( beta)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none
! parameters
  real(rtype), intent(in) :: beta
! local variables
  integer     :: i, j, ie
  real(rtype) :: two3, ui, uj, vi, vj, uav, vav, mav, kav
  real(rtype) :: txn, tyn, ds2, rds, duxa, duya, dvxa, dvya, dtxa, dtya
  real(rtype) :: duds, dvds, dtds, dudt, dvdt, dtdt
  real(rtype) :: duxf, duyf, dvxf, dvyf, dtxf, dtyf
  real(rtype) :: tauxx, tauxy, tauyy, phix, phiy
  real(rtype) :: fv(3)
! *****************************************************************************
  two3 = 2.D0/3.D0
! interior edges --------------------------------------------------------------
  do ie=1,nedint
    i = edge(1,ie)
    j = edge(2,ie)
! - average of flow variables
    ui  = cv(2,i)/cv(1,i)
    uj  = cv(2,j)/cv(1,j)
    vi  = cv(3,i)/cv(1,i)
    vj  = cv(3,j)/cv(1,j)
    uav = 0.5D0*(ui+uj)
    vav = 0.5D0*(vi+vj)
    mav = 0.5D0*(dv(6,i)+dv(6,j))
    kav = 0.5D0*(dv(7,i)+dv(7,j))
! - tangential vector (normalized)
    txn = x(j) - x(i)
    tyn = y(j) - y(i)
    ds2 = txn*txn + tyn*tyn
    rds = 1.D0/Sqrt(ds2)
    txn = txn*rds
    tyn = tyn*rds
! - average of gradients
    duxa = 0.5D0*(gradx(2,i)+gradx(2,j))
    duya = 0.5D0*(grady(2,i)+grady(2,j))
    dvxa = 0.5D0*(gradx(3,i)+gradx(3,j))
    dvya = 0.5D0*(grady(3,i)+grady(3,j))
    dtxa = 0.5D0*(gradx(5,i)+gradx(5,j))
    dtya = 0.5D0*(grady(5,i)+grady(5,j))
! - divided difference
    duds = rds*(uj-ui)
    dvds = rds*(vj-vi)
    dtds = rds*(dv(2,j)-dv(2,i))
! - tangential component - divided difference
    dudt = duxa*txn + duya*tyn - duds
    dvdt = dvxa*txn + dvya*tyn - dvds
    dtdt = dtxa*txn + dtya*tyn - dtds
! - face gradients (Eq. (5.73))
    duxf = duxa - dudt*txn
    duyf = duya - dudt*tyn
    dvxf = dvxa - dvdt*txn
    dvyf = dvya - dvdt*tyn
    dtxf = dtxa - dtdt*txn
    dtyf = dtya - dtdt*tyn
! - viscous fluxes
    tauxx = two3*mav*(2.D0*duxf-dvyf)
    tauyy = two3*mav*(2.D0*dvyf-duxf)
    tauxy =      mav*(     duyf+dvxf)
    phix  = uav*tauxx + vav*tauxy + kav*dtxf
    phiy  = uav*tauxy + vav*tauyy + kav*dtyf
    fv(1) = sij(1,ie)*tauxx + sij(2,ie)*tauxy
    fv(2) = sij(1,ie)*tauxy + sij(2,ie)*tauyy
    fv(3) = sij(1,ie)*phix  + sij(2,ie)*phiy
! - edge contributions to dissipation
    diss(2,i) = diss(2,i) + fv(1)*beta
    diss(3,i) = diss(3,i) + fv(2)*beta
    diss(4,i) = diss(4,i) + fv(3)*beta
    diss(2,j) = diss(2,j) - fv(1)*beta
    diss(3,j) = diss(3,j) - fv(2)*beta
    diss(4,j) = diss(4,j) - fv(3)*beta
  enddo
! edges to dummy nodes --------------------------------------------------------
  do ie=nedint+1,nedges
    i = edge(1,ie)
    j = edge(2,ie)
! - average of variables
    uav = 0.5D0*(cv(2,i)/cv(1,i)+cv(2,j)/cv(1,j))
    vav = 0.5D0*(cv(3,i)/cv(1,i)+cv(3,j)/cv(1,j))
    mav = 0.5D0*(dv(6,i)+dv(6,j))
    kav = 0.5D0*(dv(7,i)+dv(7,j))
! - viscous fluxes
    tauxx = two3*mav*(2.D0*gradx(2,i)-grady(3,i))
    tauyy = two3*mav*(2.D0*grady(3,i)-gradx(2,i))
    tauxy =      mav*(     grady(2,i)+gradx(3,i))
    phix  = uav*tauxx + vav*tauxy + kav*gradx(5,i)
    phiy  = uav*tauxy + vav*tauyy + kav*grady(5,i)
    fv(1) = sij(1,ie)*tauxx + sij(2,ie)*tauxy
    fv(2) = sij(1,ie)*tauxy + sij(2,ie)*tauyy
    fv(3) = sij(1,ie)*phix  + sij(2,ie)*phiy
    diss(2,i) = diss(2,i) + fv(1)*beta
    diss(3,i) = diss(3,i) + fv(2)*beta
    diss(4,i) = diss(4,i) + fv(3)*beta
  enddo
end subroutine FluxViscous
subroutine Irsmoo( ncontr,rhsold,rhsit)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : Periodic
  implicit none
! parameters
  integer     :: ncontr(:)
  real(rtype) :: rhsold(:,:), rhsit(:,:)
! local variables
  integer     :: itirs, i, j, ie, ib, ibn, ibegn, iendn
  real(rtype) :: den
! *****************************************************************************
! initialize counter (no. of contributions to a node);
! store the non-smoothed residual in "rhsold"
  do i=1,nndint
    ncontr(i)   = 0
    rhsold(1,i) = rhs(1,i)
    rhsold(2,i) = rhs(2,i)
    rhsold(3,i) = rhs(3,i)
    rhsold(4,i) = rhs(4,i)
  enddo
! Jacobi iteration ------------------------------------------------------------
  do itirs=1,nitirs
! - zero out nodal contributions
    do i=1,nndint
      rhsit(1,i) = 0D0
      rhsit(2,i) = 0D0
      rhsit(3,i) = 0D0
      rhsit(4,i) = 0D0
    enddo
! - loop over edges - first iteration => set counter
    if(itirs == 1)then
      do ie=1,nedint
        i          = edge(1,ie)
        j          = edge(2,ie)
        ncontr(i)  = ncontr(i) + 1
        ncontr(j)  = ncontr(j) + 1
        rhsit(1,i) = rhsit(1,i) + rhs(1,j)
        rhsit(2,i) = rhsit(2,i) + rhs(2,j)
        rhsit(3,i) = rhsit(3,i) + rhs(3,j)
        rhsit(4,i) = rhsit(4,i) + rhs(4,j)
        rhsit(1,j) = rhsit(1,j) + rhs(1,i)
        rhsit(2,j) = rhsit(2,j) + rhs(2,i)
        rhsit(3,j) = rhsit(3,j) + rhs(3,i)
        rhsit(4,j) = rhsit(4,j) + rhs(4,i)
      enddo
! - sum up no. of contributions at periodic nodes
      ibegn = 1
      do ib=1,nsegs
        iendn = ibound(2,ib)
        if(btype(ib)>=700 .and. btype(ib)<800)then
          do ibn=ibegn,iendn
            i         = bnode(1,ibn)
            j         = bnode(2,ibn)
            ncontr(i) = ncontr(i) + ncontr(j)
            ncontr(j) = ncontr(i)
          enddo
        endif
        ibegn = iendn + 1
      enddo
! - loop over edges - without counter
    else
      do ie=1,nedint
        i          = edge(1,ie)
        j          = edge(2,ie)
        rhsit(1,i) = rhsit(1,i) + rhs(1,j)
        rhsit(2,i) = rhsit(2,i) + rhs(2,j)
        rhsit(3,i) = rhsit(3,i) + rhs(3,j)
        rhsit(4,i) = rhsit(4,i) + rhs(4,j)
        rhsit(1,j) = rhsit(1,j) + rhs(1,i)
        rhsit(2,j) = rhsit(2,j) + rhs(2,i)
        rhsit(3,j) = rhsit(3,j) + rhs(3,i)
        rhsit(4,j) = rhsit(4,j) + rhs(4,i)
      enddo
    endif    ! itirs > 1
! - periodic boundaries
    call Periodic(rhsit)
! - new smoothed residual
    do i=1,nndint
      den      = 1.D0/(1.D0+epsirs*Real(ncontr(i)))
      rhs(1,i) = (rhsit(1,i)*epsirs+rhsold(1,i))*den
      rhs(2,i) = (rhsit(2,i)*epsirs+rhsold(2,i))*den
      rhs(3,i) = (rhsit(3,i)*epsirs+rhsold(3,i))*den
      rhs(4,i) = (rhsit(4,i)*epsirs+rhsold(4,i))*den
    enddo
  enddo    ! loop over itirs
end subroutine Irsmoo
subroutine FluxRoe2
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : FluxWalls
  implicit none
! local variables
  integer     :: i, j, ie
  real(rtype) :: rx, ry, rrho, gam1, ggm1, rl, ul, vl, pl, hl, rr, &
                 ur, vr, pr, hr, qsrl, qsrr, pav
  real(rtype) :: fc(4)
! *****************************************************************************
! initialize residual by adding artificial dissipation
  do i=1,nnodes
    rhs(1,i) = -diss(1,i)
    rhs(2,i) = -diss(2,i)
    rhs(3,i) = -diss(3,i)
    rhs(4,i) = -diss(4,i)
  enddo
! average of fluxes
  do ie=1,nedges
    i  = edge(1,ie)
    j  = edge(2,ie)
    rx = 0.5D0*(x(j)-x(i))
    ry = 0.5D0*(y(j)-y(i))
! - left & right state
    rrho = 1.D0/cv(1,i)
    gam1 = dv(4,i) - 1.D0
    ggm1 = dv(4,i)/gam1
    rl   = cv(1,i)      + lim(1,i)*(gradx(1,i)*rx+grady(1,i)*ry)
    ul   = cv(2,i)*rrho + lim(2,i)*(gradx(2,i)*rx+grady(2,i)*ry)
    vl   = cv(3,i)*rrho + lim(3,i)*(gradx(3,i)*rx+grady(3,i)*ry)
    pl   = dv(1,i)      + lim(4,i)*(gradx(4,i)*rx+grady(4,i)*ry)
    hl   = ggm1*pl/rl + 0.5D0*(ul*ul+vl*vl)
    qsrl = (ul*sij(1,ie)+vl*sij(2,ie))*rl
    rrho = 1.D0/cv(1,j)
    gam1 = dv(4,j) - 1.D0
    ggm1 = dv(4,j)/gam1
    rr   = cv(1,j)      - lim(1,j)*(gradx(1,j)*rx+grady(1,j)*ry)
    ur   = cv(2,j)*rrho - lim(2,j)*(gradx(2,j)*rx+grady(2,j)*ry)
    vr   = cv(3,j)*rrho - lim(3,j)*(gradx(3,j)*rx+grady(3,j)*ry)
    pr   = dv(1,j)      - lim(4,j)*(gradx(4,j)*rx+grady(4,j)*ry)
    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
    qsrr = (ur*sij(1,ie)+vr*sij(2,ie))*rr
! - fluxes
    pav   = 0.5D0*(pl+pr)
    fc(1) = 0.5D0*(qsrl   +qsrr  )
    fc(2) = 0.5D0*(qsrl*ul+qsrr*ur) + sij(1,ie)*pav
    fc(3) = 0.5D0*(qsrl*vl+qsrr*vr) + sij(2,ie)*pav
    fc(4) = 0.5D0*(qsrl*hl+qsrr*hr)
    rhs(1,i) = rhs(1,i) + fc(1)
    rhs(2,i) = rhs(2,i) + fc(2)
    rhs(3,i) = rhs(3,i) + fc(3)
    rhs(4,i) = rhs(4,i) + fc(4)
    rhs(1,j) = rhs(1,j) - fc(1)
    rhs(2,j) = rhs(2,j) - fc(2)
    rhs(3,j) = rhs(3,j) - fc(3)
    rhs(4,j) = rhs(4,j) - fc(4)
  enddo
! treatment of solid walls
  call FluxWalls
end subroutine FluxRoe2
!> Evaluates limiter functions (Venkatakrishnan's limiter, Eq. (5.67)).
!!
!! @param umin  minimum of U_i and of min_j U_j (U = rho, u, v, p)
!! @param umax  maximum of U_i and of max_j U_j (U = rho, u, v, p)
!!
subroutine Limiter( umin,umax)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none
! parameters
  real(rtype), intent(in) :: umin(:,:), umax(:,:)
! local variables
  integer     :: i, j, ib, ibn, ie, ibegn, iendn
  real(rtype) :: limfac3, rvolref, rx, ry, voll, eps2nl, d1minl, d1maxl, &
                 d2l, ul, vl, volr, eps2nr, d1minr, d1maxr, d2r, ur, vr, &
                 limval
  real(rtype) :: eps2(4)
! *****************************************************************************
! initialize limiter functions (1.0 = no limiting)
  do i=1,nnodes
    lim(1,i) = 1.D0
    lim(2,i) = 1.D0
    lim(3,i) = 1.D0
    lim(4,i) = 1.D0
  enddo
! normalize epsilon^2 for all limited variables (rho, u, v, p)
  limfac3 = limfac*limfac*limfac
  rvolref = 1.D0/volref**1.5D0
  eps2(1) = limfac3*limref(1)*limref(1)*rvolref
  eps2(2) = limfac3*limref(2)*limref(2)*rvolref
  eps2(3) = limfac3*limref(3)*limref(3)*rvolref
  eps2(4) = limfac3*limref(4)*limref(4)*rvolref
! evaluate limiter functions --------------------------------------------------
  do ie=1,nedint
    i = edge(1,ie)
    j = edge(2,ie)
    rx   = 0.5D0*(x(j)-x(i))
    ry   = 0.5D0*(y(j)-y(i))
    voll = vol(i)**1.5D0
    volr = vol(j)**1.5D0
! - density
    eps2nl   = eps2(1)*voll
    d1minl   = umin(1,i) - cv(1,i)
    d1maxl   = umax(1,i) - cv(1,i)
    eps2nr   = eps2(1)*volr
    d1minr   = umin(1,j) - cv(1,j)
    d1maxr   = umax(1,j) - cv(1,j)
    d2l      =  gradx(1,i)*rx + grady(1,i)*ry
    d2r      = -gradx(1,j)*rx - grady(1,j)*ry
    limval   = Venkat( d2l,d1minl,d1maxl,eps2nl)
    lim(1,i) = Min(limval,lim(1,i))
    limval   = Venkat( d2r,d1minr,d1maxr,eps2nr)
    lim(1,j) = Min(limval,lim(1,j))
! - u
    ul       = cv(2,i)/cv(1,i)
    ur       = cv(2,j)/cv(1,j)
    eps2nl   = eps2(2)*voll
    d1minl   = umin(2,i) - ul
    d1maxl   = umax(2,i) - ul
    eps2nr   = eps2(2)*volr
    d1minr   = umin(2,j) - ur
    d1maxr   = umax(2,j) - ur
    d2l      =  gradx(2,i)*rx + grady(2,i)*ry
    d2r      = -gradx(2,j)*rx - grady(2,j)*ry
    limval   = Venkat( d2l,d1minl,d1maxl,eps2nl)
    lim(2,i) = Min(limval,lim(2,i))
    limval   = Venkat( d2r,d1minr,d1maxr,eps2nr)
    lim(2,j) = Min(limval,lim(2,j))
! - v
    vl       = cv(3,i)/cv(1,i)
    vr       = cv(3,j)/cv(1,j)
    eps2nl   = eps2(3)*voll
    d1minl   = umin(3,i) - vl
    d1maxl   = umax(3,i) - vl
    eps2nr   = eps2(3)*volr
    d1minr   = umin(3,j) - vr
    d1maxr   = umax(3,j) - vr
    d2l      =  gradx(3,i)*rx + grady(3,i)*ry
    d2r      = -gradx(3,j)*rx - grady(3,j)*ry
    limval   = Venkat( d2l,d1minl,d1maxl,eps2nl)
    lim(3,i) = Min(limval,lim(3,i))
    limval   = Venkat( d2r,d1minr,d1maxr,eps2nr)
    lim(3,j) = Min(limval,lim(3,j))
! - pressure
    eps2nl   = eps2(4)*voll
    d1minl   = umin(4,i) - dv(1,i)
    d1maxl   = umax(4,i) - dv(1,i)
    eps2nr   = eps2(4)*volr
    d1minr   = umin(4,j) - dv(1,j)
    d1maxr   = umax(4,j) - dv(1,j)
    d2l      =  gradx(4,i)*rx + grady(4,i)*ry
    d2r      = -gradx(4,j)*rx - grady(4,j)*ry
    limval   = Venkat( d2l,d1minl,d1maxl,eps2nl)
    lim(4,i) = Min(limval,lim(4,i))
    limval   = Venkat( d2r,d1minr,d1maxr,eps2nr)
    lim(4,j) = Min(limval,lim(4,j))
  enddo
! periodic boundaries ---------------------------------------------------------
  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if(btype(ib)>=700 .and. btype(ib)<800)then
      do ibn=ibegn,iendn
        i = bnode(1,ibn)
        j = bnode(2,ibn)
        lim(1,i) = Min(lim(1,i),lim(1,j))
        lim(2,i) = Min(lim(2,i),lim(2,j))
        lim(3,i) = Min(lim(3,i),lim(3,j))
        lim(4,i) = Min(lim(4,i),lim(4,j))
        lim(1,j) = lim(1,i)
        lim(2,j) = lim(2,i)
        lim(3,j) = lim(3,i)
        lim(4,j) = lim(4,i)
      enddo
    endif
    ibegn = iendn + 1
  enddo
! *****************************************************************************
contains
!> Evaluates the limiter function
!!
!! @param d2     change of value to be limited (gradient * distance)
!! @param d1min  min(U_i, U_j) - U_i
!! @param d1max  max(U_i, U_j) - U_i
!! @param eps2   threshold value
!! @return       limiter value
!!
  real(rtype) function Venkat( d2,d1min,d1max,eps2)
    implicit none
    real(rtype) :: d2, d1min, d1max, eps2
    real(rtype) :: num, den
    Venkat = 1.D0
    if(d2 > 1.D-12)then
      num    = (d1max*d1max+eps2)*d2 + 2.D0*d2*d2*d1max
      den    = d2*(d1max*d1max+2.D0*d2*d2+d1max*d2+eps2)
      Venkat = num/den
    else if(d2 < -1.D-12)then
      num    = (d1min*d1min+eps2)*d2 + 2.D0*d2*d2*d1min
      den    = d2*(d1min*d1min+2.D0*d2*d2+d1min*d2+eps2)
      Venkat = num/den
    endif
  end function Venkat
end subroutine Limiter
!> Computes the convective fluxes using average of fluxes at control volume's
!! faces. The left and right fluxes are computed directly from values at nodes
!! "i" and "j". This is sufficient for a 1st-order Roe scheme.
!!
subroutine FluxRoe1
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : FluxWalls
  implicit none
! local variables
  integer     :: i, j, ie
  real(rtype) :: rhl, rhr, qsl, qsr, pav
  real(rtype) :: fc(4)
! *****************************************************************************
! initialize residual by adding artificial dissipation
  do i=1,nnodes
    rhs(1,i) = -diss(1,i)
    rhs(2,i) = -diss(2,i)
    rhs(3,i) = -diss(3,i)
    rhs(4,i) = -diss(4,i)
  enddo
! average of fluxes
  do ie=1,nedges
    i = edge(1,ie)
    j = edge(2,ie)
! - left and right rho*H, V*n
    rhl = dv(1,i) + cv(4,i)
    qsl = (cv(2,i)*sij(1,ie)+cv(3,i)*sij(2,ie))/cv(1,i)
    rhr = dv(1,j) + cv(4,j)
    qsr = (cv(2,j)*sij(1,ie)+cv(3,j)*sij(2,ie))/cv(1,j)
! - fluxes
    pav   = 0.5D0*(dv(1,i)+dv(1,j))
    fc(1) = 0.5D0*(qsl*cv(1,i)+qsr*cv(1,j))
    fc(2) = 0.5D0*(qsl*cv(2,i)+qsr*cv(2,j)) + sij(1,ie)*pav
    fc(3) = 0.5D0*(qsl*cv(3,i)+qsr*cv(3,j)) + sij(2,ie)*pav
    fc(4) = 0.5D0*(qsl*rhl    +qsr*rhr   )
    rhs(1,i) = rhs(1,i) + fc(1)
    rhs(2,i) = rhs(2,i) + fc(2)
    rhs(3,i) = rhs(3,i) + fc(3)
    rhs(4,i) = rhs(4,i) + fc(4)
    rhs(1,j) = rhs(1,j) - fc(1)
    rhs(2,j) = rhs(2,j) - fc(2)
    rhs(3,j) = rhs(3,j) - fc(3)
    rhs(4,j) = rhs(4,j) - fc(4)
  enddo
! treatment of solid walls
  call FluxWalls
end subroutine FluxRoe1
!> Computes minimum and maximum values of density, u, v, and pressure for
!! all direct neighbors "j" of node "i" (min/max_j U_j in Eq. (5.61)).
!! This is used later in Limiter to evaluate the limiter functions.
!!
!! @param umin  minimum of U_i and of min_j U_j (U = rho, u, v, p)
!! @param umax  maximum of U_i and of max_j U_j (U = rho, u, v, p)
!!
subroutine LimiterInit( umin,umax)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none
! parameters
  real(rtype), intent(out) :: umin(:,:), umax(:,:)
! local variables
  integer     :: i, j, ib, ibn, ie, ibegn, iendn
  real(rtype) :: rl, ul, vl, pl, rr, ur, vr, pr
! *****************************************************************************
! initialize with values at node i
  do i=1,nndint
    umin(1,i) = cv(1,i)
    umin(2,i) = cv(2,i)/cv(1,i)
    umin(3,i) = cv(3,i)/cv(1,i)
    umin(4,i) = dv(1,i)
    umax(1,i) = umin(1,i)
    umax(2,i) = umin(2,i)
    umax(3,i) = umin(3,i)
    umax(4,i) = umin(4,i)
  enddo
! loop over interior edges
  do ie=1,nedint
    i = edge(1,ie)
    j = edge(2,ie)
! - left state
    rl = cv(1,i)
    ul = cv(2,i)/rl
    vl = cv(3,i)/rl
    pl = dv(1,i)
! - right state
    rr = cv(1,j)
    ur = cv(2,j)/rr
    vr = cv(3,j)/rr
    pr = dv(1,j)
! - neighbors of node i
    umin(1,i) = Min(umin(1,i),rr)
    umin(2,i) = Min(umin(2,i),ur)
    umin(3,i) = Min(umin(3,i),vr)
    umin(4,i) = Min(umin(4,i),pr)
    umax(1,i) = Max(umax(1,i),rr)
    umax(2,i) = Max(umax(2,i),ur)
    umax(3,i) = Max(umax(3,i),vr)
    umax(4,i) = Max(umax(4,i),pr)
! - neighbors of node j
    umin(1,j) = Min(umin(1,j),rl)
    umin(2,j) = Min(umin(2,j),ul)
    umin(3,j) = Min(umin(3,j),vl)
    umin(4,j) = Min(umin(4,j),pl)
    umax(1,j) = Max(umax(1,j),rl)
    umax(2,j) = Max(umax(2,j),ul)
    umax(3,j) = Max(umax(3,j),vl)
    umax(4,j) = Max(umax(4,j),pl)
  enddo
! periodic boundaries
  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if(btype(ib)>=700 .and. btype(ib)<800)then
      do ibn=ibegn,iendn
        i = bnode(1,ibn)
        j = bnode(2,ibn)
        umin(1,i) = Min(umin(1,i),umin(1,j))
        umin(2,i) = Min(umin(2,i),umin(2,j))
        umin(3,i) = Min(umin(3,i),umin(3,j))
        umin(4,i) = Min(umin(4,i),umin(4,j))
        umax(1,i) = Max(umax(1,i),umax(1,j))
        umax(2,i) = Max(umax(2,i),umax(2,j))
        umax(3,i) = Max(umax(3,i),umax(3,j))
        umax(4,i) = Max(umax(4,i),umax(4,j))
        umin(1,j) = umin(1,i)
        umin(2,j) = umin(2,i)
        umin(3,j) = umin(3,i)
        umin(4,j) = umin(4,i)
        umax(1,j) = umax(1,i)
        umax(2,j) = umax(2,i)
        umax(3,j) = umax(3,i)
        umax(4,j) = umax(4,i)
      enddo
    endif
    ibegn = iendn + 1
  enddo
end subroutine LimiterInit
!> Corrects face vectors of edges which represent the symmetry boundary.
!! The reason is that there must be no component in the direction normal
!! to the boundary for the fluxes to be computed correctly. Subroutine
!! also changes the boundary type to reflect the orientation of the symmetry
!! boundary: 501 for symmetry along x-direction, 502 along y-direction.
!!
!! @param marker  temporary memory for a node marker
!!
subroutine FaceVectorsSymm( marker)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  implicit none
! parameters
  integer :: marker(:)
! local variables
  integer     :: i, j, ib, ibf, ibn, ie, ibegf, ibegn, iendf, iendn
  real(rtype) :: sx, sy
! *****************************************************************************
! mark nodes at symmetry boundaries
  do i=1,nndint
    marker(i) = -1
  enddo
  ibegf = 1
  ibegn = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    iendn = ibound(2,ib)
    if(btype(ib)>=500 .and. btype(ib)<600)then
      sx = 0D0
      sy = 0D0
      do ibf=ibegf,iendf
        sx = sx + sbf(1,ibf)
        sy = sy + sbf(2,ibf)
      enddo
      if(Abs(sx) > Abs(sy))then
        btype(ib) = 501      ! symmetry at x=const. plane
      else
        btype(ib) = 502      ! symmetry at y=const. plane
      endif
      do ibn=ibegn,iendn
        marker(bnode(1,ibn)) = btype(ib) - 500    ! store symmetry axis
      enddo
    endif
    ibegf = iendf + 1
    ibegn = iendn + 1
  enddo
! correct face vectors
  do ie=1,nedint
    i = edge(1,ie)
    j = edge(2,ie)
    if(marker(i)/=-1 .and. marker(j)/=-1)then
      if(marker(i) < 2)then     ! x=const. plane
        sij(1,ie) = 0D0
      else                        ! y=const. plane
        sij(2,ie) = 0D0
      endif
    endif
  enddo
end subroutine FaceVectorsSymm
!> Computes reference values of limited variables (density, u, v, pressure)
!! and of the control volumes. The reference values are used to normalize
!! variables within the limiter functions (Roe's upwind scheme).
!!
subroutine LimiterRefvals
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none
! local variables
  real(rtype) :: gam1, rgas, temp, rho, cs, mach
! *****************************************************************************
! reference volume (= largest control volume)
  volref = Maxval( vol)
! reference density, velocity and pressure
  gam1 = gamma - 1.D0
  rgas = gam1*cpgas/gamma
! external flow
  if(kflow == "E")then
    limref(1) = rhoinf
    limref(2) = Sqrt(uinf*uinf+vinf*vinf)
    limref(3) = limref(2)
    limref(4) = pinf
! internal flow
  else
    temp      = ttinl*(pout/ptinl)**(gam1/gamma)
    rho       = pout/(rgas*temp)
    cs        = Sqrt(gamma*pout/rho)
    mach      = Sqrt(2.D0*((ttinl/temp)-1.D0)/gam1)
    limref(1) = rho
    limref(2) = mach*cs
    limref(3) = limref(2)
    limref(4) = pout
  endif
end subroutine LimiterRefvals
!> Prints out error message and stops program execution.
!!
!! @param message error message
!!
subroutine ErrorMessage( message)
  implicit none
! parameters
  character(*), intent(in) :: message
! *****************************************************************************
  write(*,"(A,A)")" Error: ",Trim( message)
  stop
end subroutine ErrorMessage
!> Generates temporary lists with nodes of an edge (niedge, iedge).
!! Computes total number of edges (interior + dummy). The edge lists
!! are used in the subroutines EdgesFinalize and InitMetrics.
!!
!! @param niedge  pointer from a node to iedge()
!! @param iedge   linked list of edge endpoints:
!!                @li (1,*) = point j of edge (i,j)
!!                @li (2,*) = next point j which is also connected to i;
!!                            if <0 - no further connections
!!                @li (3,*) = pointer to edge() - used in InitMetrics to associate
!!                            face vector sij() with the correct edge
!!
subroutine EdgesInitialize( niedge,iedge)
  use ModGeometry
  use ModNumerics
  use ModInterfaces, only : ErrorMessage
  implicit none
! parameters
  integer, intent(out) :: niedge(:), iedge(:,:)
! local variables
  integer :: d, cedge, cedge2, mxedges
  integer :: i, j, ic, ie, n
! *****************************************************************************
  mxedges = ubound(iedge,2)  ! max. possible number of edges
! reset all pointers
  do i=1,nndint
    niedge(i) = -777
  enddo
  do ie=1,mxedges
    iedge(1,ie) = -777
    iedge(2,ie) = -777
    iedge(3,ie) = -777
  enddo
! loop over nodes of all triangles
  nedint = 0
  do n=1,3
! - loop over triangles
    do ic=1,ntria
      i = tria(n,ic)
      if(n < 3)then
        j = tria(n+1,ic)
      else
        j = tria(1,ic)
      endif
      if(i > j)then  ! lower index first
        d = i
        i = j
        j = d
      endif
      if(niedge(i) < 0)then
! ----- define new edge
        nedint = nedint + 1
        if(nedint > mxedges)then
          call ErrorMessage("max. number of edges reached")
        endif
        niedge(i)       = nedint
        iedge(1,nedint) = j
      else
! ----- insert node "j" into list of adjacent nodes
        cedge = niedge(i)
10      continue
          if(iedge(1,cedge) == j) goto 20
          cedge2 = iedge(2,cedge)
          if(cedge2 < 0)then
            nedint = nedint + 1
            if(nedint > mxedges)then
              call ErrorMessage("max. number of edges reached")
            endif
            iedge(2,cedge) = nedint
            iedge(1,nedint) = j
            goto 20
          endif
          cedge = cedge2
        goto 10
20      continue
      endif
    enddo ! loop over triangles
  enddo   ! loop over nodes of triangles
! set total no. of edges (add edges to dummy nodes)
  nedges = nedint + (nnodes-nndint)
end subroutine EdgesInitialize
!> Computes mass flow at inlet and mass flow ratio between inlet and outlet
!! boundaries. Contributions are summed up by looping over ALL inlet and
!! outlet boundaries.
!!
subroutine Massflow
  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModPlotQuant
  implicit none
! local variables
  integer     :: ibegf, iendf, n1, n2, in, out
  integer     :: ib, ibf
  real(rtype) :: sx, sy, massin, massout, mass
! *****************************************************************************
! initialize mass flow before summing up the contributions
  massin  = 0D0
  massout = 0D0
  mflow   = 0D0
  mfratio = 0D0
  in  = 0  ! flow into the domain (=1)
  out = 0  ! flow out of the domain (=1)
! loop over boundaries searching for inlet / outlet
  ibegf = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    if(btype(ib)>=100 .and. btype(ib)<300)then
      do ibf=ibegf,iendf
        n1   = bface(1,ibf)
        n2   = bface(2,ibf)
        sx   = sbf(1,ibf)
        sy   = sbf(2,ibf)
        mass = 0.5D0*((cv(2,n1)+cv(2,n2))*sx+(cv(3,n1)+cv(3,n2))*sy)
        if(btype(ib)>=100 .and. btype(ib)<200)then
! ------- inflow
          massin  = massin - mass
          in = 1
        else
! ------- outflow
          massout = massout + mass
          out = 1
        endif
      enddo
    endif
    ibegf = iendf + 1
  enddo
! mass flow and ratio
  if(in == 1)then
    mflow   = massin
    mfratio = massout/massin
  endif
  if(in==0 .and. out==1)then
    mflow   = massout
    mfratio = 1.D0
  endif
end subroutine Massflow
!> Generates final edge list (interior edges, edges to dummy nodes) using
!! the temporary lists "niedge" and "iedge".
!!
!! @param niedge  pointer from a node to iedge()
!! @param iedge   linked list of edge endpoints:
!!                @li (1,*) = point j of edge (i,j)
!!                @li (2,*) = next point j which is also connected to i;
!!                            if <0 - no further connections
!!                @li (3,*) = pointer to edge() - used in InitMetrics to associate
!!                            face vector sij() with the correct edge
!!
subroutine EdgesFinalize( niedge,iedge)
  use ModGeometry
  use ModNumerics
  use ModInterfaces, only : ErrorMessage
  implicit none
! parameters
  integer :: niedge(:), iedge(:,:)
! local variables
  integer :: errFlag, i, ibn, ie, cedge
! *****************************************************************************
  allocate( edge(2,nedges),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for edge()")
  ie = 0  ! edge counter
  do i=1,nndint
! - loop over all grid nodes
    cedge = niedge(i)
    if(cedge > 0)then
10    continue
! ----- loop over all edges adjacent to node "i"
        ie             = ie + 1
        edge(1,ie)     = i
        edge(2,ie)     = iedge(1,cedge)
        iedge(3,cedge) = ie                 ! we need it in InitMetrics
        cedge          = iedge(2,cedge)     ! next adjacent edge
        if(cedge < 0) goto 20
      goto 10
20    continue
    endif
  enddo
  if(ie /= nedint)then
    call ErrorMessage("did not get the correct number of interior edges")
  endif
! add edges to dummy nodes;
! store 'dummy' edges in "bnode(3,*)"
  do ibn=1,nbnodes
    if(bnode(3,ibn) == -1)then      ! dummy node here (see DummyNodes)
      ie           = ie + 1
      edge(1,ie)   = bnode(1,ibn)     ! boundary node first
      edge(2,ie)   = bnode(2,ibn)     ! dummy node second
      bnode(3,ibn) = ie
    endif
  enddo
  if(ie /= nedges)then
    call ErrorMessage("did not get the correct number of dummy edges")
  endif
end subroutine EdgesFinalize
!> Stores indexes of dummy nodes in "bnode" -> index of boundary node,
!! index of dummy node. Adds number of dummy nodes to number of physical
!! (interior) nodes.
!!
subroutine DummyNodes
  use ModGeometry
  use ModInterfaces, only : ErrorMessage
  implicit none
! local variables
  character(chrlen) :: msg
  logical :: flag
  integer :: errFlag, ibegf, iendf, ibegn, iendn, itype
  integer :: i, ib, ibf, idn
  integer, allocatable :: marker(:)  ! node marker
! *****************************************************************************
  allocate( marker(nndint),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate temporary marker")
! loop over boundary segments
  ibegf = 1
  ibegn = 1
  idn   = 0  ! counter of dummy nodes
  do ib=1,nsegs
    iendf = ibound(1,ib)
    iendn = ibound(2,ib)
    itype = btype(ib)
    flag  = .false.  ! true for inlet/outlet/far-field
    if(itype>=100 .and. itype<200) flag = .true.
    if(itype>=200 .and. itype<300) flag = .true.
    if(itype>=600 .and. itype<700) flag = .true.
    if(itype<700 .or. itype>=800)then   ! NOT a periodic boundary
! --- reset node marker
      do i=1,nndint
        marker(i) = -777
      enddo
! --- loop over faces of boundary "ib" and mark nodes
      do ibf=ibegf,iendf
        marker(bface(1,ibf)) = 1
        marker(bface(2,ibf)) = 1
      enddo
! --- store node indexes in "bnode";
!     count dummy nodes (idn) for inlet/outlet/far-field boundary
      do i=1,nndint
        if(marker(i) == 1)then           ! must be on boundary
          if(ibegn > nbnodes)then        ! check dimension
            call ErrorMessage("max. no. of boundary nodes exceeded")
          endif
          if(flag)then                   ! *** inlet/outlet/far-field
            idn            = idn + 1
            bnode(1,ibegn) = i             ! index of boundary node
            bnode(2,ibegn) = nndint + idn  ! index of dummy node
            bnode(3,ibegn) = -1            ! set in "EdgesFinalize"
          else                             ! *** other boundary type
            bnode(1,ibegn) = i             ! index of boundary node
          endif
          ibegn = ibegn + 1                ! count boundary nodes
        endif
      enddo
! --- check no. of boundary nodes
      if(ibegn-1 .ne. iendn)then
        write(msg,1000) itype,iendn,ibegn-1
        call ErrorMessage(msg)
      endif
    endif  ! periodic?
! - reset pointers to faces and to nodes
    ibegf = iendf + 1
    ibegn = iendn + 1
  enddo  ! ib
  deallocate( marker)
! set total number of nodes (add dummy nodes)
  nnodes = nndint + idn
1000 format("no. of nodes for boundary ",I3," is wrong. It should be ",I5, &
            " but it is ",I5)
end subroutine DummyNodes
!> Computes upwind dissipation according to 2nd-order Roe's flux-difference
!! splitting scheme. Values are extrapolated to the dual faces using Venkat's
!! limiter function. The dissipation terms are multiplied by a preconditioning
!! matrix for low Mach numbers (see Eqs. (9.54), (9.82)).
!!
!! @param beta  coefficient for mixing new and old dissipation values
!!
subroutine DissipRoe2Prec( beta)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : CompTheta, RightEigenvec, MatprodTp1_P1, MatVecProd5
  implicit none
! parameters
  real(rtype), intent(in) :: beta
! local variables
  integer     :: i, j, ie
  real(rtype) :: beta5, ds, rx, ry, gam1, ggm1, rgas, rrho, rl, ul, vl, pl, &
                 tl, hl, rr, ur, vr, pr, tr, hr, rav, dd, dd1, uav, vav, hav, &
                 pav, tav, cav, q2a, uv, h1, h2, h5, delta, eabs1, eabs2, &
                 eabs5, a1, a4, a5, cc2, ra1g, rhop, rhoT, hT, theta
  real(rtype) :: fd(5), nVec(3), gtpMat(5,5), tp1Mat(5,5)
  real(rtype) :: wVec(5), wpVec(5), wrlVec(5), dumVec(5)
! *****************************************************************************
  beta5 = 0.5D0*beta
  nVec(3)   = 0D0   ! no 3rd dimension here
  wrlVec(4) = 0D0
  wVec(4)   = 0D0
  wpVec(4)  = 0D0
  do ie=1,nedges
    i       = edge(1,ie)
    j       = edge(2,ie)
    ds      = Sqrt(sij(1,ie)*sij(1,ie)+sij(2,ie)*sij(2,ie))
    nVec(1) = sij(1,ie)/ds
    nVec(2) = sij(2,ie)/ds
    rx      = 0.5D0*(x(j)-x(i))
    ry      = 0.5D0*(y(j)-y(i))
! - left & right state
    rrho = 1.D0/cv(1,i)
    gam1 = dv(4,i) - 1.D0
    ggm1 = dv(4,i)/gam1
    rgas = gam1*dv(5,i)/dv(4,i)
    rl   = cv(1,i)      + lim(1,i)*(gradx(1,i)*rx+grady(1,i)*ry)
    ul   = cv(2,i)*rrho + lim(2,i)*(gradx(2,i)*rx+grady(2,i)*ry)
    vl   = cv(3,i)*rrho + lim(3,i)*(gradx(3,i)*rx+grady(3,i)*ry)
    pl   = dv(1,i)      + lim(4,i)*(gradx(4,i)*rx+grady(4,i)*ry)
    tl   = pl/(rgas*rl)
    hl   = ggm1*pl/rl + 0.5D0*(ul*ul+vl*vl)
    rrho = 1.D0/cv(1,j)
    gam1 = dv(4,j) - 1.D0
    ggm1 = dv(4,j)/gam1
    rgas = gam1*dv(5,j)/dv(4,j)
    rr   = cv(1,j)      - lim(1,j)*(gradx(1,j)*rx+grady(1,j)*ry)
    ur   = cv(2,j)*rrho - lim(2,j)*(gradx(2,j)*rx+grady(2,j)*ry)
    vr   = cv(3,j)*rrho - lim(3,j)*(gradx(3,j)*rx+grady(3,j)*ry)
    pr   = dv(1,j)      - lim(4,j)*(gradx(4,j)*rx+grady(4,j)*ry)
    tr   = pr/(rgas*rr)
    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
! - Roe's average
    rav  = Sqrt(rl*rr)
    gam1 = 0.5D0*(dv(4,i)+dv(4,j)) - 1.D0
    dd   = rav/rl
    dd1  = 1.D0/(1.D0+dd)
    uav  = (ul+dd*ur)*dd1
    vav  = (vl+dd*vr)*dd1
    pav  = (pl+dd*pr)*dd1
    tav  = (tl+dd*tr)*dd1
    hav  = (hl+dd*hr)*dd1
    q2a  = uav*uav + vav*vav
    cav  = Sqrt(gam1*(hav-0.5D0*q2a))
    uv   = uav*nVec(1) + vav*nVec(2)
! - preconditioning
    wrlVec(1) = rr - rl
    wrlVec(2) = rr*ur - rl*ul
    wrlVec(3) = rr*vr - rl*vl
    wrlVec(5) = (rr*hr-pr) - (rl*hl-pl)
    wVec(1) = rav
    wVec(2) = rav*uav
    wVec(3) = rav*vav
    wVec(5) = rav*hav - pav
    wpVec(1) = pav
    wpVec(2) = uav
    wpVec(3) = vav
    wpVec(5) = tav
    rhop  =  rav/pav
    rhoT  = -rav/tav
    hT    = 0.5D0*(dv(5,i)+dv(5,j))
    theta = CompTheta( gam1+1.D0,cav,q2a)
    delta = epsentr*cav
    h2    = Abs(uv)
    eabs2 = EntropyCorr2Prec( h2,delta)
    a1    = wVec(1)*rhop*hT + rhoT
    ra1g  = 1.D0/(wVec(1)*theta*hT + rhoT)
    a4    = a1*ra1g
    a5    = wVec(1)*hT*ra1g
    cc2   = 0.5D0*Sqrt((uv*uv)*((a4-1.D0)*(a4-1.D0))+4.D0*a5)
    h2    = 0.5D0*(a4+1.D0)*uv
    h1    = Abs(h2 + cc2)
    h5    = Abs(h2 - cc2)
    eabs1 = EntropyCorr2Prec( h1,delta)
    eabs5 = EntropyCorr2Prec( h5,delta)
    call RightEigenvec(wVec,wpVec,nVec,uv,hav,theta,rhop,rhoT,0D0,hT,gtpMat)
    call MatprodTp1_P1(wVec,wpVec,nVec,uv,hav,theta,rhop,rhoT,0D0,hT,q2a,tp1Mat)
! - upwind fluxes
    call MatVecProd5(tp1Mat,wrlVec,dumVec)
    dumVec(1) = dumVec(1)*eabs2
    dumVec(2) = dumVec(2)*eabs2
    dumVec(3) = dumVec(3)*eabs2
    dumVec(4) = dumVec(4)*eabs1
    dumVec(5) = dumVec(5)*eabs5
    call MatVecProd5(gtpMat,dumVec,fd)
! - edge contributions to dissipation
    ds        = ds*beta5
    diss(1,i) = diss(1,i) + fd(1)*ds
    diss(2,i) = diss(2,i) + fd(2)*ds
    diss(3,i) = diss(3,i) + fd(3)*ds
    diss(4,i) = diss(4,i) + fd(5)*ds
    diss(1,j) = diss(1,j) - fd(1)*ds
    diss(2,j) = diss(2,j) - fd(2)*ds
    diss(3,j) = diss(3,j) - fd(3)*ds
    diss(4,j) = diss(4,j) - fd(5)*ds
  enddo
! *****************************************************************************
contains
!> Evaluates entropy correction function
!!
!! @param z  value to be corrected
!! @param d  threshold value of the correction
!! @return   corrected value of z
!!
  real(rtype) function EntropyCorr2Prec( z,d)
    implicit none
    real(rtype) :: z, d
    if(z > d)then
      EntropyCorr2Prec = z
    else
      EntropyCorr2Prec = 0.5D0*(z*z+d*d)/d
    endif
  end function EntropyCorr2Prec
end subroutine DissipRoe2Prec
!> Computes upwind dissipation according to 2nd-order Roe's flux-difference
!! splitting scheme. Values are extrapolated to the dual faces using Venkat's
!! limiter function.
!!
!! @param beta  coefficient for mixing new and old dissipation values
!!
subroutine DissipRoe2( beta)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none
! parameters
  real(rtype), intent(in) :: beta
! local variables
  integer     :: i, j, ie
  real(rtype) :: beta5, ds, nx, ny, rx, ry, gam1, ggm1, rrho, rl, ul, vl, &
                 pl, hl, rr, ur, vr, pr, hr, rav, dd, dd1, uav, vav, hav, &
                 q2a, c2a, cav, uv, du, h1, h2, h3, h4, h5, delta, eabs1, &
                 eabs2, eabs4
  real(rtype) :: fd(4)
! *****************************************************************************
  beta5 = 0.5D0*beta
  do ie=1,nedges
    i  = edge(1,ie)
    j  = edge(2,ie)
    ds = Sqrt(sij(1,ie)*sij(1,ie)+sij(2,ie)*sij(2,ie))
    nx = sij(1,ie)/ds
    ny = sij(2,ie)/ds
    rx = 0.5D0*(x(j)-x(i))
    ry = 0.5D0*(y(j)-y(i))
! - left & right state
    rrho = 1.D0/cv(1,i)
    gam1 = dv(4,i) - 1.D0
    ggm1 = dv(4,i)/gam1
    rl   = cv(1,i)      + lim(1,i)*(gradx(1,i)*rx+grady(1,i)*ry)
    ul   = cv(2,i)*rrho + lim(2,i)*(gradx(2,i)*rx+grady(2,i)*ry)
    vl   = cv(3,i)*rrho + lim(3,i)*(gradx(3,i)*rx+grady(3,i)*ry)
    pl   = dv(1,i)      + lim(4,i)*(gradx(4,i)*rx+grady(4,i)*ry)
    hl   = ggm1*pl/rl + 0.5D0*(ul*ul+vl*vl)
    rrho = 1.D0/cv(1,j)
    gam1 = dv(4,j) - 1.D0
    ggm1 = dv(4,j)/gam1
    rr   = cv(1,j)      - lim(1,j)*(gradx(1,j)*rx+grady(1,j)*ry)
    ur   = cv(2,j)*rrho - lim(2,j)*(gradx(2,j)*rx+grady(2,j)*ry)
    vr   = cv(3,j)*rrho - lim(3,j)*(gradx(3,j)*rx+grady(3,j)*ry)
    pr   = dv(1,j)      - lim(4,j)*(gradx(4,j)*rx+grady(4,j)*ry)
    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
! - Roe's average
    rav   = Sqrt(rl*rr)
    gam1  = 0.5D0*(dv(4,i)+dv(4,j)) - 1.D0
    dd    = rav/rl
    dd1   = 1.D0/(1.D0+dd)
    uav   = (ul+dd*ur)*dd1
    vav   = (vl+dd*vr)*dd1
    hav   = (hl+dd*hr)*dd1
    q2a   = 0.5D0*(uav*uav+vav*vav)
    c2a   = gam1*(hav-q2a)
    cav   = Sqrt(c2a)
    uv    = uav*nx + vav*ny
    du    = (ur-ul)*nx + (vr-vl)*ny
! - eigenvalues
    h1    = Abs(uv - cav)
    h2    = Abs(uv)
    h4    = Abs(uv + cav)
    delta = epsentr*h4
    eabs1 = EntropyCorr2( h1,delta)
    eabs2 = EntropyCorr2( h2,delta)
    eabs4 = EntropyCorr2( h4,delta)
! - upwind fluxes
    h1 = rav*cav*du
    h2 = eabs1*(pr-pl - h1)/(2.D0*c2a)
    h3 = eabs2*(rr-rl - (pr-pl)/c2a)
    h4 = eabs2*rav
    h5 = eabs4*(pr-pl + h1)/(2.D0*c2a)
    fd(1) = h2 + h3 + h5
    fd(2) = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-du*nx) + h5*(uav+cav*nx)
    fd(3) = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-du*ny) + h5*(vav+cav*ny)
    fd(4) = h2*(hav-cav*uv) + h3*q2a + h4*(uav*(ur-ul)+vav*(vr-vl)-uv*du) + &
            h5*(hav+cav*uv)
! - edge contributions to dissipation
    ds        = ds*beta5
    diss(1,i) = diss(1,i) + fd(1)*ds
    diss(2,i) = diss(2,i) + fd(2)*ds
    diss(3,i) = diss(3,i) + fd(3)*ds
    diss(4,i) = diss(4,i) + fd(4)*ds
    diss(1,j) = diss(1,j) - fd(1)*ds
    diss(2,j) = diss(2,j) - fd(2)*ds
    diss(3,j) = diss(3,j) - fd(3)*ds
    diss(4,j) = diss(4,j) - fd(4)*ds
  enddo
! *****************************************************************************
contains
!> Evaluates entropy correction function
!!
!! @param z  value to be corrected
!! @param d  threshold value of the correction
!! @return   corrected value of z
!!
  real(rtype) function EntropyCorr2( z,d)
    implicit none
    real(rtype) :: z, d
    if(z > d)then
      EntropyCorr2 = z
    else
      EntropyCorr2 = 0.5D0*(z*z+d*d)/d
    endif
  end function EntropyCorr2
end subroutine DissipRoe2
!> Computes upwind dissipation according to 1st-order Roe's flux-difference
!! splitting scheme. The dissipation terms are multiplied by a preconditioning
!! matrix for low Mach numbers (see Eqs. (9.54), (9.82)).
!!
!! @param beta  coefficient for mixing new and old dissipation values
!!
subroutine DissipRoe1Prec( beta)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : CompTheta, RightEigenvec, MatprodTp1_P1, MatVecProd5
  implicit none
! parameters
  real(rtype), intent(in) :: beta
! local variables
  integer     :: i, j, ie
  real(rtype) :: beta5, ds, gam1, rrho, rl, ul, vl, pl, tl, hl, &
                 rr, ur, vr, pr, tr, hr, rav, dd, dd1, uav, vav, pav, tav, &
                 hav, cav, q2a, uv, h1, h2, h5, delta, eabs1, eabs2, eabs5, &
                 a1, a4, a5, cc2, ra1g, rhop, rhoT, hT, theta
  real(rtype) :: fd(5), nVec(3), gtpMat(5,5), tp1Mat(5,5)
  real(rtype) :: wVec(5), wpVec(5), wrlVec(5), dumVec(5)
! *****************************************************************************
  beta5 = 0.5D0*beta
  nVec(3)   = 0D0   ! no 3rd dimension here
  wrlVec(4) = 0D0
  wVec(4)   = 0D0
  wpVec(4)  = 0D0
  do ie=1,nedges
    i       = edge(1,ie)
    j       = edge(2,ie)
    ds      = Sqrt(sij(1,ie)*sij(1,ie)+sij(2,ie)*sij(2,ie))
    nVec(1) = sij(1,ie)/ds
    nVec(2) = sij(2,ie)/ds
! - left & right state
    rrho = 1.D0/cv(1,i)
    rl   = cv(1,i)
    ul   = cv(2,i)*rrho
    vl   = cv(3,i)*rrho
    pl   = dv(1,i)
    tl   = dv(2,i)
    hl   = (pl+cv(4,i))*rrho
    rrho = 1.D0/cv(1,j)
    rr   = cv(1,j)
    ur   = cv(2,j)*rrho
    vr   = cv(3,j)*rrho
    pr   = dv(1,j)
    tr   = dv(2,j)
    hr   = (pr+cv(4,j))*rrho
! - Roe's average
    rav  = Sqrt(rl*rr)
    gam1 = 0.5D0*(dv(4,i)+dv(4,j)) - 1.D0
    dd   = rav/rl
    dd1  = 1.D0/(1.D0+dd)
    uav  = (ul+dd*ur)*dd1
    vav  = (vl+dd*vr)*dd1
    pav  = (pl+dd*pr)*dd1
    tav  = (tl+dd*tr)*dd1
    hav  = (hl+dd*hr)*dd1
    q2a  = uav*uav + vav*vav
    cav  = Sqrt(gam1*(hav-0.5D0*q2a))
    uv   = uav*nVec(1) + vav*nVec(2)
! - preconditioning
    wrlVec(1) = rr - rl
    wrlVec(2) = cv(2,j) - cv(2,i)
    wrlVec(3) = cv(3,j) - cv(3,i)
    wrlVec(5) = cv(4,j) - cv(4,i)
    wVec(1) = rav
    wVec(2) = rav*uav
    wVec(3) = rav*vav
    wVec(5) = rav*hav - pav
    wpVec(1) = pav
    wpVec(2) = uav
    wpVec(3) = vav
    wpVec(5) = tav
    rhop  =  rav/pav
    rhoT  = -rav/tav
    hT    = 0.5D0*(dv(5,i)+dv(5,j))
    theta = CompTheta( gam1+1.D0,cav,q2a)
    delta = epsentr*cav
    h2    = Abs(uv)
    eabs2 = EntropyCorr1Prec( h2,delta)
    a1    = wVec(1)*rhop*hT + rhoT
    ra1g  = 1.D0/(wVec(1)*theta*hT + rhoT)
    a4    = a1*ra1g
    a5    = wVec(1)*hT*ra1g
    cc2   = 0.5D0*Sqrt((uv*uv)*((a4-1.D0)*(a4-1.D0))+4.D0*a5)
    h2    = 0.5D0*(a4+1.D0)*uv
    h1    = Abs(h2 + cc2)
    h5    = Abs(h2 - cc2)
    eabs1 = EntropyCorr1Prec( h1,delta)
    eabs5 = EntropyCorr1Prec( h5,delta)
    call RightEigenvec(wVec,wpVec,nVec,uv,hav,theta,rhop,rhoT,0D0,hT,gtpMat)
    call MatprodTp1_P1(wVec,wpVec,nVec,uv,hav,theta,rhop,rhoT,0D0,hT,q2a,tp1Mat)
! - upwind fluxes
    call MatVecProd5(tp1Mat,wrlVec,dumVec)
    dumVec(1) = dumVec(1)*eabs2
    dumVec(2) = dumVec(2)*eabs2
    dumVec(3) = dumVec(3)*eabs2
    dumVec(4) = dumVec(4)*eabs1
    dumVec(5) = dumVec(5)*eabs5
    call MatVecProd5(gtpMat,dumVec,fd)
! - edge contributions to dissipation
    ds        = ds*beta5
    diss(1,i) = diss(1,i) + fd(1)*ds
    diss(2,i) = diss(2,i) + fd(2)*ds
    diss(3,i) = diss(3,i) + fd(3)*ds
    diss(4,i) = diss(4,i) + fd(5)*ds
    diss(1,j) = diss(1,j) - fd(1)*ds
    diss(2,j) = diss(2,j) - fd(2)*ds
    diss(3,j) = diss(3,j) - fd(3)*ds
    diss(4,j) = diss(4,j) - fd(5)*ds
  enddo
! *****************************************************************************
contains
!> Evaluates entropy correction function
!!
!! @param z  value to be corrected
!! @param d  threshold value of the correction
!! @return   corrected value of z
!!
  real(rtype) function EntropyCorr1Prec( z,d)
    implicit none
    real(rtype) :: z, d
    if(z > d)then
      EntropyCorr1Prec = z
    else
      EntropyCorr1Prec = 0.5D0*(z*z+d*d)/d
    endif
  end function EntropyCorr1Prec
end subroutine DissipRoe1Prec
!> Computes upwind dissipation according to 1st-order Roe's flux-difference
!! splitting scheme.
!!
!! @param beta  coefficient for mixing new and old dissipation values
!!
subroutine DissipRoe1( beta)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none
! parameters
  real(rtype), intent(in) :: beta
! local variables
  integer     :: i, j, ie
  real(rtype) :: beta5, ds, nx, ny, gam1, rrho, rl, ul, vl, pl, hl, &
                 rr, ur, vr, pr, hr, rav, dd, dd1, uav, vav, hav, q2a, &
                 c2a, cav, uv, du, h1, h2, h3, h4, h5, delta, eabs1, &
                 eabs2, eabs4
  real(rtype) :: fd(4)
! *****************************************************************************
  beta5 = 0.5D0*beta
  do ie=1,nedges
    i  = edge(1,ie)
    j  = edge(2,ie)
    ds = Sqrt(sij(1,ie)*sij(1,ie)+sij(2,ie)*sij(2,ie))
    nx = sij(1,ie)/ds
    ny = sij(2,ie)/ds
! - left & right state
    rrho = 1.D0/cv(1,i)
    rl   = cv(1,i)
    ul   = cv(2,i)*rrho
    vl   = cv(3,i)*rrho
    pl   = dv(1,i)
    hl   = (pl+cv(4,i))*rrho
    rrho = 1.D0/cv(1,j)
    rr   = cv(1,j)
    ur   = cv(2,j)*rrho
    vr   = cv(3,j)*rrho
    pr   = dv(1,j)
    hr   = (pr+cv(4,j))*rrho
! - Roe's average
    rav   = Sqrt(rl*rr)
    gam1  = 0.5D0*(dv(4,i)+dv(4,j)) - 1.D0
    dd    = rav/rl
    dd1   = 1.D0/(1.D0+dd)
    uav   = (ul+dd*ur)*dd1
    vav   = (vl+dd*vr)*dd1
    hav   = (hl+dd*hr)*dd1
    q2a   = 0.5D0*(uav*uav+vav*vav)
    c2a   = gam1*(hav-q2a)
    cav   = Sqrt(c2a)
    uv    = uav*nx + vav*ny
    du    = (ur-ul)*nx + (vr-vl)*ny
! - eigenvalues
    h1    = Abs(uv - cav)
    h2    = Abs(uv)
    h4    = Abs(uv + cav)
    delta = epsentr*h4
    eabs1 = EntropyCorr1( h1,delta)
    eabs2 = EntropyCorr1( h2,delta)
    eabs4 = EntropyCorr1( h4,delta)
! - upwind fluxes
    h1 = rav*cav*du
    h2 = eabs1*(pr-pl - h1)/(2.D0*c2a)
    h3 = eabs2*(rr-rl - (pr-pl)/c2a)
    h4 = eabs2*rav
    h5 = eabs4*(pr-pl + h1)/(2.D0*c2a)
    fd(1) = h2 + h3 + h5
    fd(2) = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-du*nx) + h5*(uav+cav*nx)
    fd(3) = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-du*ny) + h5*(vav+cav*ny)
    fd(4) = h2*(hav-cav*uv) + h3*q2a + h4*(uav*(ur-ul)+vav*(vr-vl)-uv*du) + &
            h5*(hav+cav*uv)
! - edge contributions to dissipation
    ds        = ds*beta5
    diss(1,i) = diss(1,i) + fd(1)*ds
    diss(2,i) = diss(2,i) + fd(2)*ds
    diss(3,i) = diss(3,i) + fd(3)*ds
    diss(4,i) = diss(4,i) + fd(4)*ds
    diss(1,j) = diss(1,j) - fd(1)*ds
    diss(2,j) = diss(2,j) - fd(2)*ds
    diss(3,j) = diss(3,j) - fd(3)*ds
    diss(4,j) = diss(4,j) - fd(4)*ds
  enddo
! *****************************************************************************
contains
!> Evaluates entropy correction function
!!
!! @param z  value to be corrected
!! @param d  threshold value of the correction
!! @return   corrected value of z
!!
  real(rtype) function EntropyCorr1( z,d)
    implicit none
    real(rtype) :: z, d
    if(z > d)then
      EntropyCorr1 = z
    else
      EntropyCorr1 = 0.5D0*(z*z+d*d)/d
    endif
  end function EntropyCorr1
end subroutine DissipRoe1
!> Computes values of dependent variables (pressure, temperature, speed
!! of sound, specific heat ratio, specific heat coeff. at const. pressure)
!! from conservative variables at all grid points. Additionally, laminar
!! viscosity and heat conductivity coefficients are computed in the case
!! of viscous flow.
!!
subroutine DependentVarsAll
  use ModDataTypes
  use ModGeometry
  use ModPhysics
  implicit none
! local variables
  integer     :: i
  real(rtype) :: gam1, rgas, g1cp, rhoq, s1, s2, s12, rat, cppr
! *****************************************************************************
  gam1 = gamma - 1.D0
  rgas = gam1*cpgas/gamma
  g1cp = gam1*cpgas
! Euler equations
  if(kequs == "E")then
    do i=1,nnodes
      rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
      dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
      dv(2,i) = dv(1,i)/(rgas*cv(1,i))
      dv(3,i) = Sqrt(g1cp*dv(2,i))
      dv(4,i) = gamma
      dv(5,i) = cpgas
    enddo
! Navier-Stokes equations
  else
    s1   = 110D0
    s2   = 288.16D0
    s12  = 1.D0 + s1/s2
    cppr = cpgas/prlam
    do i=1,nnodes
      rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
      dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
      dv(2,i) = dv(1,i)/(rgas*cv(1,i))
      dv(3,i) = Sqrt(g1cp*dv(2,i))
      dv(4,i) = gamma
      dv(5,i) = cpgas
      rat     = Sqrt(dv(2,i)/s2)*s12/(1.D0+s1/dv(2,i))
      dv(6,i) = refvisc*rat
      dv(7,i) = dv(6,i)*cppr
    enddo
  endif
end subroutine DependentVarsAll
!> Computes values of dependent variables (pressure, temperature, speed
!! of sound, specific heat ratio, specific heat coeff. at const. pressure)
!! from conservative variables at the node i. Additionally, laminar
!! viscosity and heat conductivity coefficients are computed in the case
!! of viscous flow.
!!
!! @param i  node index
!!
subroutine DependentVarsOne( i)
  use ModDataTypes
  use ModPhysics
  implicit none
! parameters
  integer, intent(in) :: i
! local variables
  real(rtype) :: gam1, rgas, g1cp, rhoq, s1, s2, s12, rat
! *****************************************************************************
  gam1 = gamma - 1.D0
  rgas = gam1*cpgas/gamma
  g1cp = gam1*cpgas
! Euler equations
  if(kequs == "E")then
    rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
    dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
    dv(2,i) = dv(1,i)/(rgas*cv(1,i))
    dv(3,i) = Sqrt(g1cp*dv(2,i))
    dv(4,i) = gamma
    dv(5,i) = cpgas
! Navier-Stokes equations
  else
    s1      = 110D0
    s2      = 288.16D0
    s12     = 1.D0 + s1/s2
    rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
    dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
    dv(2,i) = dv(1,i)/(rgas*cv(1,i))
    dv(3,i) = Sqrt(g1cp*dv(2,i))
    dv(4,i) = gamma
    dv(5,i) = cpgas
    rat     = Sqrt(dv(2,i)/s2)*s12/(1.D0+s1/dv(2,i))
    dv(6,i) = refvisc*rat
    dv(7,i) = dv(6,i)*(cpgas/prlam)
  endif
end subroutine DependentVarsOne
!> Monitors the convergence, prints it out and stores it in a file. For external
!! flow, it also prints out the lift, the drag and the moment coefficients. For
!! internal flow, it prints out the mass flow and the mass flow ratio.
!!
subroutine Convergence
  use ModDataTypes
  use ModFiles
  use ModControl
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : Forces, Massflow
  implicit none
! local variables
  integer     :: i, idr
  real(rtype) :: dr, drmax
! *****************************************************************************
! compute the residual
  drho  = 0D0
  drmax = 0D0
  do i=1,nndint
    dr   = cv(1,i) - cvold(1,i)
    drho = drho + dr*dr
    if(Abs(dr) >= drmax)then
      drmax = Abs(dr)
      idr   = i
    endif
  enddo
  if(iter == 1)then
    drho1 = Sqrt(drho) + 1.D-32
    drho  = 1.D0
  else
    drho  = Sqrt(drho)/drho1
  endif
! compute forces & moments (external flow)
  if(kflow == "E")then
    call Forces
! compute mass flow and mass flow ratio (internal flow)
  else
    call Massflow
  endif
! print out / store
  if(kflow == "E")then
!!  write(ifConv,1000) iter,Log10(drho),drmax,idr,cl,cd,cm
    write(*,1005) iter,Log10(drho),drmax,idr,cl,cd,cm
  else
!!  write(ifConv,1010) iter,Log10(drho),drmax,idr,mflow,mfratio
    write(*,1015) iter,Log10(drho),drmax,idr,mflow,mfratio
  endif
1000  format(I5,1P,2X,E12.5,2X,E12.5,0P,I8,1P,3(2X,E12.5))
1005  format(I5,1P,2X,E11.4,2X,E10.4,0P,I8,1P,3(2X,E10.3))
1010  format(I5,1P,2X,E12.5,2X,E12.5,0P,I8,1X,1P,2(2X,E12.5))
1015  format(I5,1P,2X,E11.4,2X,E10.4,0P,I8,1X,1P,2(2X,E10.4))
end subroutine Convergence
!> Checks metrics by computing min. and max. volume, and by computing
!! the sum of face vectors for each control volume (should be zero).
!!
!! @param work  work space (used to store the sum of face vectors
!!              of a control volume - fvecSum)
!!
subroutine CheckMetrics( work)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  implicit none
! parameters
  real(rtype) :: work(:)
! local variables
  integer     :: i, j, ib, ie, ibegf, iendf, ibegn, iendn
  real(rtype) :: volmin, volmax, s, smax
  real(rtype), allocatable :: fvecSum(:,:)
! *****************************************************************************
!!!allocate( work(2*(nconv*nnodes)),stat=errFlag)
!!!fvecSum = Reshape( work,(/2, nndint/))
  allocate(fvecSum(2,nndint))
! obtain mix. and max. control volume
  volmin = Minval( vol)
  volmax = Maxval( vol)
! compute sum of face vectors for each control volume (fvecSum)
  do i=1,nndint
    fvecSum(1,i) = 0D0
    fvecSum(2,i) = 0D0
  enddo
  do ie=1,nedint
    i            = edge(1,ie)
    j            = edge(2,ie)
    fvecSum(1,i) = fvecSum(1,i) + sij(1,ie)
    fvecSum(2,i) = fvecSum(2,i) + sij(2,ie)
    fvecSum(1,j) = fvecSum(1,j) - sij(1,ie)
    fvecSum(2,j) = fvecSum(2,j) - sij(2,ie)
  enddo
  ibegf = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    if(btype(ib)<700 .or. btype(ib)>=800)then   ! boundary faces (non-periodic)
      do ie=ibegf,iendf
        i            = bface(1,ie)
        j            = bface(2,ie)
        fvecSum(1,i) = fvecSum(1,i) + 0.5D0*sbf(1,ie)
        fvecSum(2,i) = fvecSum(2,i) + 0.5D0*sbf(2,ie)
        fvecSum(1,j) = fvecSum(1,j) + 0.5D0*sbf(1,ie)
        fvecSum(2,j) = fvecSum(2,j) + 0.5D0*sbf(2,ie)
      enddo
    endif
    ibegf = iendf + 1
  enddo
  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if(btype(ib)>=700 .and. btype(ib)<800)then   ! periodic nodes
      do ie=ibegn,iendn
        i            = bnode(1,ie)
        j            = bnode(2,ie)
        fvecSum(1,i) = fvecSum(1,i) + fvecSum(1,j)
        fvecSum(2,i) = fvecSum(2,i) + fvecSum(2,j)
        fvecSum(1,j) = fvecSum(1,i)
        fvecSum(2,j) = fvecSum(2,i)
      enddo
    endif
    ibegn = iendn + 1
  enddo
! compute maximum of the sum of face vectors
  smax = -1.D+32
  do i=1,nndint
    s    = Sqrt(fvecSum(1,i)**2+fvecSum(2,i)**2)
    smax = Max(smax,s)
  enddo
! print info
  write(*,1000) smax,volmin,volmax
1000  format(" max. sum(S) = ",E11.4,/," min. volume = ",E11.4,/, &
             " max. volume = ",E11.4)
end subroutine CheckMetrics
!> Adds both parts of variable var(:,:) at all periodic boundaries.
!!
!! @param var  variable to update
!!
subroutine Periodic( var)
  use ModDataTypes
  use ModGeometry
  implicit none
! parameters
  real(rtype) :: var(:,:)
! local variables
  integer :: i, j, ib, ibn, ibegn, iendn, n
! *****************************************************************************
  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if(btype(ib)>=700 .and. btype(ib)<800)then
      do n=1,Ubound(var,1)
        do ibn=ibegn,iendn
          i        = bnode(1,ibn)
          j        = bnode(2,ibn)
          var(n,i) = var(n,i) + var(n,j)
          var(n,j) = var(n,i)
        enddo
      enddo
    endif
    ibegn = iendn + 1
  enddo
end subroutine Periodic
!> Sets boundary conditions at dummy points. In a first loop, b.c.'s at
!! inlet, outlet and far-field are specified. In a second loop, b.c.'s are
!! set for all solid walls.
!!
!! @param work  work space for temporary variables (used by far-field b.c.)
!!
subroutine BoundaryConditions( work)
  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModInterfaces, only : BcondFarfield, BcondInflow, BcondOutflow, &
                            BcondWallns, ErrorMessage
  implicit none
! parameters
  real(rtype) :: work(:)
! local variables
  integer :: ib, ibegn, iendn, itype, wdim
! *****************************************************************************
! loop over all boundaries with the exception of walls
  ibegn = 1
  do ib=1,nsegs
    itype = btype(ib)
    iendn = ibound(2,ib)
! - inflow
    if(itype>=100 .and. itype<200)then
      call BcondInflow(ibegn,iendn)
! - outflow
    else if(itype>=200 .and. itype<300)then
      call BcondOutflow(ibegn,iendn)
! - far-field
    else if(itype>=600 .and. itype<700)then
      wdim = Ubound(work,1)
      if((4*nbnodes) > wdim)then
        call ErrorMessage("insufficient work space in BoundaryConditions")
      endif
      call BcondFarfield(ibegn,iendn, &
                          work( 1           :  nbnodes), &
                          work((1+  nbnodes):2*nbnodes), &
                          work((1+2*nbnodes):3*nbnodes), &
                          work((1+3*nbnodes):4*nbnodes))
    endif
    ibegn = iendn + 1
  enddo ! ib
! solid walls (treated last because they should dominate)
  ibegn = 1
  do ib=1,nsegs
    itype = btype(ib)
    iendn = ibound(2,ib)
! - viscous (no-slip) wall - if Navier-Stokes equations solved
    if(itype>=300 .and. itype<400 .and. kequs=="N")then
      call BcondWallns(ibegn,iendn)
    endif
    ibegn = iendn + 1
  enddo ! ib
end subroutine BoundaryConditions
!> Writes out selected quantities for the whole flow field in Vis2D format.
!!
subroutine PlotFlow
  use ModDataTypes
  use ModControl
  use ModFiles
  use ModGeometry
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : ErrorMessage
  implicit none
! local variables
  character(chrlen) :: fname
  integer     :: errFlag, nquant, i, m,i1,i2,i3
  real(rtype) :: rrho, u, v, e, press, temp, c, mach, ttot, ptot, machis, &
                 ptloss, pratio, ptotinf, gam1, ggm1, visc,xc,yc,uc,vc
  real(rtype) :: varout(mxqfield+2),u1,u2,u3,v1,v2,v3,x1,x2,x3,y1,y2,y3
! *****************************************************************************
  ptotinf = 0D0
  if(kflow == "E")then
    gam1    = gamma - 1.D0
    ggm1    = gamma/gam1
    ptotinf = pinf*(1.D0+0.5D0*gam1*machinf*machinf)**ggm1
  endif
  open(1,file='unstruct2d.v2d',status='unknown')
  open(2,file='unstruct2d.2dv',status='unknown')
  open(3,file='unstruct2d.plt',status='unknown')
  write(3,3)
3 format('TITLE="2D Fluid Flow"')
  write(3,4)
4 format('VARIABLES="X" "Y" "P" "U" "V"')
  write(3,5)nndint,ntria
5 format('ZONE T="pressure" N=',I5,', E=',I5,', F=FEPOINT, ET=TRIANGLE')
  write(3,6)
6 format('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)')
  write(2,'(I5)')nndint
  do i=1,nndint
!   rrho =1D0/cv(1,i)
!   u    =cv(2,i)*rrho
!   v    =cv(3,i)*rrho
!   e    =cv(4,i)*rrho
    press=dv(1,i)
!   temp =dv(2,i)
!   c    =dv(3,i)
!   gam1 =dv(4,i)-1D0
!   ggm1 =dv(4,i)/gam1
!   if(kequs.eq."N")then
!     visc=dv(6,i)
!   else
!     visc=0D0
!   endif
!   mach=Sqrt(u*u+v*v)/c
!   ttot=(e+press*rrho)/dv(5,i)
!   ptot=press*(ttot/temp)**ggm1
!   if(kflow.eq."E")then
!     ptloss=1D0-ptot/ptotinf
!     pratio=ptotinf/press
!   else
!     ptloss=1D0-ptot/ptinl
!     pratio=ptinl/press
!   endif
!   machis=(pratio**(1D0/ggm1)-1D0)*2D0/gam1
!   machis=Max(machis,0D0)
!   machis=Sqrt(machis)
    write(2,7)x(i),y(i),(press-pinf)/1D3,0D0,0D0
    write(3,7)x(i),y(i),(press-pinf)/1D3,0D0,0D0
7   format(5E13.5)
  enddo
  write(2,'(I5)')ntria
  do i=1,ntria
    write(2,'(3I6)')tria(1,i),tria(2,i),tria(3,i)
    write(3,'(3I6)')tria(1,i),tria(2,i),tria(3,i)
  enddo
  write(3,8)
8 format('ZONE T="velocities", F=POINT')
  do i=1,ntria
    i1=tria(1,i)
    x1=x(i1)
    y1=y(i1)
    rrho=1D0/cv(1,i1)
    u1   =cv(2,i1)*rrho
    v1   =cv(3,i1)*rrho
    i2=tria(2,i)
    x2=x(i2)
    y2=y(i2)
    rrho=1D0/cv(1,i2)
    u2   =cv(2,i2)*rrho
    v2   =cv(3,i2)*rrho
    i3=tria(3,i)
    x3=x(i3)
    y3=y(i3)
    rrho=1D0/cv(1,i3)
    u3   =cv(2,i3)*rrho
    v3   =cv(3,i3)*rrho
    xc=(x1+x2+x3)/3D0
    yc=(y1+y2+y3)/3D0
    uc=(u1+u2+u3)/3D0
    vc=(v1+v2+v3)/3D0
    write(1,7)xc,yc,    uc,vc
    write(3,7)xc,yc,0D0,uc,vc
  enddo
  close(3)
  close(2)
  close(1)
end subroutine PlotFlow
!> Applies no-slip wall boundary condition. Adiabatic walls
!! only are assumed (velocity components are zeroed out).
!!
!! @param ibegn  indirect pointer to first node of the boundary
!! @param iendn  indirect pointer to last node of the boundary
!!
subroutine BcondWallns( ibegn,iendn)
  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModInterfaces, only : DependentVarsOne
  implicit none
! parameters
  integer, intent(in) :: ibegn, iendn
! local variables
  integer :: ib, ibn
! *****************************************************************************
  do ib=ibegn,iendn
    ibn       = bnode(1,ib)   ! boundary node
    cv(2,ibn) = 0D0
    cv(3,ibn) = 0D0
    call DependentVarsOne(ibn)
  enddo
end subroutine BcondWallns
!> Writes out selected values at walls and symmetry boundaries in Vis2D format.
!!
subroutine PlotSurfaces
  use ModDataTypes
  use ModControl
  use ModFiles
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : ErrorMessage
  implicit none
! local variables
  character(chrlen) :: fname
  integer     :: errFlag, itype, ibegf, iendf, ibegn, iendn, ibf1, ibf2, &
                 nquant, nsurfs
  integer     :: i, ib, ibf, ibn, m
  real(rtype) :: rrho, u, v, e, press, temp, c, ptot, ttot, mach, machis, &
                 ptloss, pratio, ptotinf, gam1, ggm1
  real(rtype) :: cf, cp, visc, sx, sy, ds, sxn, syn, grdnx, grdny, grdnn, &
                 dvdnx, dvdny, dvdna, sgn
  real(rtype) :: varout(mxquant+2)
! *****************************************************************************
  ptotinf = 0D0
  if(kflow == "E")then
    gam1    = gamma - 1.D0
    ggm1    = gamma/gam1
    ptotinf = pinf*(1.D0+0.5D0*gam1*machinf*machinf)**ggm1
  endif
! open plot file
!!write(fname,"(A,I5.5,A)") Trim(fnSurf),iter,".v2d"
!!open(unit=ifSurf, file=fname, status="unknown", action="write", iostat=errFlag)
!!if (errFlag /= 0) call ErrorMessage("cannot open plot file(surfaces)")
! header
  nquant = 0
  do m=1,mxquant
    if(lquant(m) == "Y") nquant = nquant + 1
  enddo
  nsurfs = 0  ! no. of surfaces to output
  do ib=1,nsegs
    itype = btype(ib)
    if((itype>=300 .and. itype<600) .or. &
        (itype>=800 .and. itype<900)) nsurfs = nsurfs + 1
  enddo
!!write(ifSurf,1000) Trim(title),nsurfs,nquant+2
! names of variables
!!do m=1,mxquant
!!  if(lquant(m) == "Y")then
!!    write(ifSurf,"(A)") Trim(cquant(m))
!!  endif
!!enddo
! compute quantities & write'em out
  ibegf = 1
  ibegn = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    iendn = ibound(2,ib)
    itype = btype(ib)
    if((itype>=300 .and. itype<600) .or. &
        (itype>=800 .and. itype<900))then
!!    write(ifSurf,1010) iendn-ibegn+1,Trim(bname(ib))
      do ibn=ibegn,iendn
        i = bnode(1,ibn)
        varout(1) = x(i)
        varout(2) = y(i)
        rrho  = 1.D0/cv(1,i)
        u     = cv(2,i)*rrho
        v     = cv(3,i)*rrho
        e     = cv(4,i)*rrho
        press = dv(1,i)
        temp  = dv(2,i)
        c     = dv(3,i)
        gam1  = dv(4,i) - 1.D0
        ggm1  = dv(4,i)/gam1
        if(kequs == "N")then
          ibf1 = -1
          ibf2 = -1
          do ibf=ibegf,iendf
            if(bface(1,ibf) == i)then
              if(ibf1 < 0)then
                ibf1 = ibf
              else
                ibf2 = ibf
              endif
            endif
            if(bface(2,ibf) == i)then
              if(ibf1 < 0)then
                ibf1 = ibf
              else
                ibf2 = ibf
              endif
            endif
          enddo
          if(ibf2 < 0) ibf2 = ibf1
          visc  = dv(6,i)
          sx    = -0.5D0*(sbf(1,ibf1)+sbf(1,ibf2))   ! to point inside
          sy    = -0.5D0*(sbf(2,ibf1)+sbf(2,ibf2))
          ds    = Sqrt(sx*sx+sy*sy)
          sxn   = sx/ds
          syn   = sy/ds
          grdnx = gradx(2,i)*sxn + grady(2,i)*syn
          grdny = gradx(3,i)*sxn + grady(3,i)*syn
          grdnn = grdnx*sxn + grdny*syn
          dvdnx = grdnx - grdnn*sxn
          dvdny = grdny - grdnn*syn
          if(grdnx > grdny)then    ! to get somehow the main flow
            sgn = Sign(1.D0,grdnx)
          else
            sgn = Sign(1.D0,grdny)
          endif
          dvdna = Sqrt(dvdnx*dvdnx+dvdny*dvdny)
          if(kflow == "E")then
            cf = 2.D0*sgn*visc*dvdna/(rhoinf*qinf*qinf)
          else
            cf = 2.D0*sgn*visc*dvdna/(refrho*refvel*refvel)
          endif
        else
          visc = 0D0
          cf   = 0D0
        endif
        mach  = Sqrt(u*u+v*v)/c
        ttot  = (e+press*rrho)/dv(5,i)
        ptot  = press*(ttot/temp)**ggm1
        if(kflow == "E")then
          ptloss = 1.D0 - ptot/ptotinf
          pratio = ptotinf/press
          cp     = 2.*(pinf-press)/(rhoinf*qinf*qinf)
        else
          ptloss = 1.D0 - ptot/ptinl
          pratio = ptinl/press
          cp     = 2.*((p12rat*pout)-press)/(refrho*refvel*refvel)
        endif
        machis = (pratio**(1.D0/ggm1)-1.D0)*2.D0/gam1
        machis = Max(machis, 0D0)
        machis = Sqrt(machis)
! ----- store quantities in varout()
        nquant = 2
        do m=1,mxquant
          if(lquant(m) == "Y")then
            nquant = nquant + 1
! --------- density
            if(m == 1)then
              varout(nquant) = cv(1,i)
! --------- u-velocity
            else if(m == 2)then
              varout(nquant) = u
! --------- v-velocity
            else if(m == 3)then
              varout(nquant) = v
! --------- static pressure
            else if(m == 4)then
              varout(nquant) = press
! --------- total pressure
            else if(m == 5)then
              varout(nquant) = ptot
! --------- static temperature
            else if(m == 6)then
              varout(nquant) = temp
! --------- total temperature
            else if(m == 7)then
              varout(nquant) = ttot
! --------- local Mach number
            else if(m == 8)then
              varout(nquant) = mach
! --------- isentropic Mach number
            else if(m == 9)then
              varout(nquant) = machis
! --------- total pressure loss
            else if(m == 10)then
              varout(nquant) = ptloss
! --------- laminar viscosity coefficient
            else if(m == 11)then
              varout(nquant) = visc
! --------- skin friction coefficient
            else if(m == 12)then
              varout(nquant) = cf
! --------- pressure coefficient
            else if(m == 13)then
              varout(nquant) = cp
            endif
          endif
        enddo
!!      write(ifSurf,1020) (varout(m), m=1,nquant)
      enddo ! node
    endif ! itype
    ibegf = iendf + 1
    ibegn = iendn + 1
  enddo ! ib
!!close(unit=ifSurf)
1000  format(A,/,"1",/,"Boundaries",/,I3,I3,/,"x [m]",/,"y [m]")
1010  format(I6," 0",/,"0 0 0",/,A)
1020  format(1P,20E16.8)
end subroutine PlotSurfaces
!> Applies outflow boundary condition to dummy points. Sub- or supersonic
!! outflow is possible.
!!
!! @param ibegn  indirect pointer to first node of the boundary
!! @param iendn  indirect pointer to last node of the boundary
!!
subroutine BcondOutflow( ibegn,iendn)
  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModInterfaces, only : DependentVarsOne
  implicit none
! parameters
  integer, intent(in) :: ibegn, iendn
! local variables
  integer     :: ib, ibn, idn, ie
  real(rtype) :: gam1, ds, sxn, syn, rrho, u, v, rrhoc, deltp, &
                 rhob, ub, vb, vnd, q, mach
! *****************************************************************************
  do ib=ibegn,iendn
    ibn = bnode(1,ib)         ! boundary node
    idn = bnode(2,ib)         ! dummy node
    ie  = bnode(3,ib)         ! edge to dummy node
    ds  = Sqrt(sij(1,ie)**2+sij(2,ie)**2)
    sxn = sij(1,ie)/ds
    syn = sij(2,ie)/ds
    gam1  = dv(4,ibn) - 1.D0
    rrho  = 1.D0/cv(1,ibn)
    u     = cv(2,ibn)*rrho
    v     = cv(3,ibn)*rrho
    q     = Sqrt(u*u+v*v)
    mach  = q/dv(3,ibn)
    if(mach < 1.D0)then
      rrhoc = rrho/dv(3,ibn)
      deltp = dv(1,ibn) - pout
      rhob  = cv(1,ibn) - deltp/(dv(3,ibn)*dv(3,ibn))
      ub    = u + sxn*deltp*rrhoc
      vb    = v + syn*deltp*rrhoc
! --- special treatment to prevent "deltp" from changing the sign
!     of velocity components. This may happen for very small u, v.
      vnd = ub*sxn + vb*syn
      if(vnd < 0D0)then
        ub = Sign(1.D0,u)*Max(Abs(ub),Abs(u))
        vb = Sign(1.D0,v)*Max(Abs(vb),Abs(v))
      endif
      cv(1,idn) = rhob
      cv(2,idn) = rhob*ub
      cv(3,idn) = rhob*vb
      cv(4,idn) = pout/gam1 + 0.5D0*rhob*(ub*ub+vb*vb)
    else
      cv(1,idn) = cv(1,ibn)
      cv(2,idn) = cv(2,ibn)
      cv(3,idn) = cv(3,ibn)
      cv(4,idn) = cv(4,ibn)
    endif
    call DependentVarsOne(idn)
  enddo  ! ib
end subroutine BcondOutflow
!> Computes the preconditioning parameter theta.
!!
!! @param gam  ratio of specific heats
!! @param c    speed of sound
!! @param q2   total velocity squared
!!
function CompTheta( gam,c,q2)
  use ModDataTypes
  use ModNumerics
  use ModPhysics
  implicit none
! parameters
  real(rtype), intent(in) :: gam, c, q2
! result
  real(rtype) :: CompTheta
! local variables
  real(rtype) :: mach2, mref2, beta
! *****************************************************************************
  mach2 = q2/(c*c)
  mref2 = Max(Min(mach2,1.D0),precoeff*machinf*machinf)
  beta  = mref2/(1.D0+(gam-1.D0)*mref2)
  CompTheta = 1.D0/(beta*c*c)
end function CompTheta
! =============================================================================
!> Computes transformation matrix from primitive to conservative
!! variables P (equivalent to the preconditioning matrix G).
!!
!! @param wvec   vector of conservative variables
!! @param wpvec  vector of primitive variables
!! @param H      total enthalpy
!! @param theta  preconditioning parameter
!! @param rhoT   derivative of density wrp. to temperature
!! @param hp     derivative of enthalpy wrp. to pressure
!! @param hT     derivative of enthalpy wrp. to temperature
!! @param pmat   matrix P
!!
subroutine Prim2Cons( wvec,wpvec,H,theta,rhoT,hp,hT,pmat)
  use ModDataTypes
  implicit none
! parameters
  real(rtype), intent( in) :: H, theta, rhoT, hp, hT
  real(rtype), intent( in) :: wvec(5), wpvec(5)
  real(rtype), intent(out) :: pmat(5,5)
! *****************************************************************************
  pmat(1,1) = theta
  pmat(2,1) = theta*wpvec(2)
  pmat(3,1) = theta*wpvec(3)
  pmat(4,1) = theta*wpvec(4)
  pmat(5,1) = theta*H - 1.D0 - wvec(1)*hp
  pmat(1,2) = 0D0
  pmat(2,2) = wvec(1)
  pmat(3,2) = 0D0
  pmat(4,2) = 0D0
  pmat(5,2) = wvec(2)
  pmat(1,3) = 0D0
  pmat(2,3) = 0D0
  pmat(3,3) = wvec(1)
  pmat(4,3) = 0D0
  pmat(5,3) = wvec(3)
  pmat(1,4) = 0D0
  pmat(2,4) = 0D0
  pmat(3,4) = 0D0
  pmat(4,4) = wvec(1)
  pmat(5,4) = wvec(4)
  pmat(1,5) = rhoT
  pmat(2,5) = rhoT*wpvec(2)
  pmat(3,5) = rhoT*wpvec(3)
  pmat(4,5) = rhoT*wpvec(4)
  pmat(5,5) = rhoT*H + wvec(1)*hT
end subroutine Prim2Cons
! =============================================================================
!> Computes transformation matrix from conservative to primitive variables
!! P^-1 (equivalent to the inverse of the preconditioning matrix G^-1).
!!
!! @param wvec   vector of conservative variables
!! @param wpvec  vector of primitive variables
!! @param H      total enthalpy
!! @param q2     total velocity squared
!! @param theta  preconditioning parameter
!! @param rhoT   derivative of density wrp. to temperature
!! @param hp     derivative of enthalpy wrp. to pressure
!! @param hT     derivative of enthalpy wrp. to temperature
!! @param pmat1  inverse of matrix P
!!
subroutine Cons2Prim( wvec,wpvec,H,q2,theta,rhoT,hp,hT,pmat1)
  use ModDataTypes
  implicit none
! parameters
  real(rtype), intent( in) :: H, q2, theta, rhoT, hp, hT
  real(rtype), intent( in) :: wvec(5), wpvec(5)
  real(rtype), intent(out) :: pmat1(5,5)
! local variables
  real(rtype) :: rrho, Hq2, ra1, d1, d2
! *****************************************************************************
  rrho = 1.D0/wvec(1)
  Hq2  = H - q2
  ra1  = 1.D0/(wvec(1)*theta*hT + rhoT*(1.D0-wvec(1)*hp))
  d1   = rhoT*ra1
  d2   = theta*ra1
  pmat1(1,1) = ra1*(wvec(1)*hT+rhoT*Hq2)
  pmat1(2,1) = -wpvec(2)*rrho
  pmat1(3,1) = -wpvec(3)*rrho
  pmat1(4,1) = -wpvec(4)*rrho
  pmat1(5,1) = ra1*(1.D0-theta*Hq2-wvec(1)*hp)
  pmat1(1,2) = d1*wpvec(2)
  pmat1(2,2) = rrho
  pmat1(3,2) = 0D0
  pmat1(4,2) = 0D0
  pmat1(5,2) = -d2*wpvec(2)
  pmat1(1,3) = d1*wpvec(3)
  pmat1(2,3) = 0D0
  pmat1(3,3) = rrho
  pmat1(4,3) = 0D0
  pmat1(5,3) = -d2*wpvec(3)
  pmat1(1,4) = d1*wpvec(4)
  pmat1(2,4) = 0D0
  pmat1(3,4) = 0D0
  pmat1(4,4) = rrho
  pmat1(5,4) = -d2*wpvec(4)
  pmat1(1,5) = -d1
  pmat1(2,5) = 0D0
  pmat1(3,5) = 0D0
  pmat1(4,5) = 0D0
  pmat1(5,5) = d2
end subroutine Cons2Prim
! =============================================================================
!> Computes matrix of left eigenvectors of G^-1*A_c,p (i.e., (T_p)^-1).
!!
!! @param wvec   vector of conservative variables
!! @param wpvec  vector of primitive variables
!! @param nvec   components of the unit normal vector
!! @param V      contravariant velocity
!! @param theta  preconditioning parameter
!! @param rhop   derivative of density wrp. to pressure
!! @param rhoT   derivative of density wrp. to temperature
!! @param hp     derivative of enthalpy wrp. to pressure
!! @param hT     derivative of enthalpy wrp. to temperature
!! @param evl    matrix of left eigenvectors (T_p)^-1
!!
subroutine LeftEigenvec( wvec,wpvec,nvec,V,theta,rhop,rhoT,hp,hT,evl)
  use ModDataTypes
  implicit none
! parameters
  real(rtype), intent( in) :: V, theta, rhop, rhoT, hp, hT
  real(rtype), intent( in) :: wvec(5), wpvec(5), nvec(3)
  real(rtype), intent(out) :: evl(5,5)
! local variables
  real(rtype) :: a1, a1g, ra1g, a4, a5, a6, a7, cc, h1, h2, h3
! *****************************************************************************
  a1   = wvec(1)*rhop*hT + rhoT*(1.D0-wvec(1)*hp)
  a1g  = wvec(1)*theta*hT + rhoT*(1.D0-wvec(1)*hp)
  ra1g = 1.D0/a1g
  a4   = a1*ra1g
  a5   = wvec(1)*hT*ra1g
  cc   = 0.5D0*Sqrt((V*V)*(a4-1.D0)*(a4-1.D0)+4.D0*a5)
  a6   = wvec(1)*a5
  a7   = (V*(a4-1.D0))/(4.D0*cc)
  h1   = rhoT*ra1g
  h2   = a6/cc
  h3   = a6/wpvec(5)
  evl(1,1) = -h1*nvec(1)
  evl(2,1) = -h1*nvec(2)
  evl(3,1) = -h1*nvec(3)
  evl(4,1) =  0.5D0 + a7
  evl(5,1) =  0.5D0 - a7
  evl(1,2) =  0D0
  evl(2,2) = -h2*nvec(3)
  evl(3,2) =  h2*nvec(2)
  evl(4,2) =  0.5D0*h2*nvec(1)
  evl(5,2) = -0.5D0*h2*nvec(1)
  evl(1,3) =  h2*nvec(3)
  evl(2,3) =  0D0
  evl(3,3) = -h2*nvec(1)
  evl(4,3) =  0.5D0*h2*nvec(2)
  evl(5,3) = -0.5D0*h2*nvec(2)
  evl(1,4) = -h2*nvec(2)
  evl(2,4) =  h2*nvec(1)
  evl(3,4) =  0D0
  evl(4,4) =  0.5D0*h2*nvec(3)
  evl(5,4) = -0.5D0*h2*nvec(3)
  evl(1,5) = -h3*nvec(1)
  evl(2,5) = -h3*nvec(2)
  evl(3,5) = -h3*nvec(3)
  evl(4,5) =  0D0
  evl(5,5) =  0D0
end subroutine LeftEigenvec
! =============================================================================
!> Computes matrix of right eigenvectors of G^-1*A_c,p multiplied
!! by the preconditioning matrix G (i.e., G*T_p).
!!
!! @param wvec   vector of conservative variables
!! @param wpvec  vector of primitive variables
!! @param nvec   components of the unit normal vector
!! @param V      contravariant velocity
!! @param H      total enthalpy
!! @param theta  preconditioning parameter
!! @param rhop   derivative of density wrp. to pressure
!! @param rhoT   derivative of density wrp. to temperature
!! @param hp     derivative of enthalpy wrp. to pressure
!! @param hT     derivative of enthalpy wrp. to temperature
!! @param evr    matrix of right eigenvectors multiplied by G (G*T_p)
!!
subroutine RightEigenvec( wvec,wpvec,nvec,V,H,theta,rhop,rhoT,hp,hT,evr)
  use ModDataTypes
  implicit none
! parameters
  real(rtype), intent( in) :: V, H, theta, rhop, rhoT, hp, hT
  real(rtype), intent( in) :: wvec(5), wpvec(5), nvec(3)
  real(rtype), intent(out) :: evr(5,5)
! local variables
  real(rtype) :: a1, a1g, ra1g, a4, a5, a8, a9, a10, cc, h1, h2, h3
! *****************************************************************************
  a1   = wvec(1)*rhop*hT + rhoT*(1.D0-wvec(1)*hp)
  a1g  = wvec(1)*theta*hT + rhoT*(1.D0-wvec(1)*hp)
  ra1g = 1.D0/a1g
  a4   = a1*ra1g
  a5   = wvec(1)*hT*ra1g
  cc   = 0.5D0*Sqrt((V*V)*(a4-1.D0)*(a4-1.D0)+4.D0*a5)
  a8   = rhoT*wpvec(5)/wvec(1)
  a9   = -0.5D0*V*(a4-1.D0)
  a10  = a8*H + hT*wpvec(5)
  h1   = a9 + cc
  h2   = a9 - cc
  h3   = a1g/(wvec(1)*hT)
  evr(1,1) = -a8*nvec(1)
  evr(2,1) = -a8*wpvec(2)*nvec(1)
  evr(3,1) =  cc*nvec(3) - a8*wpvec(3)*nvec(1)
  evr(4,1) = -cc*nvec(2) - a8*wpvec(4)*nvec(1)
  evr(5,1) =  cc*(wpvec(3)*nvec(3)-wpvec(4)*nvec(2)) - a10*nvec(1)
  evr(1,2) = -a8*nvec(2)
  evr(2,2) = -cc*nvec(3) - a8*wpvec(2)*nvec(2)
  evr(3,2) = -a8*wpvec(3)*nvec(2)
  evr(4,2) =  cc*nvec(1) - a8*wpvec(4)*nvec(2)
  evr(5,2) =  cc*(wpvec(4)*nvec(1)-wpvec(2)*nvec(3)) - a10*nvec(2)
  evr(1,3) = -a8*nvec(3)
  evr(2,3) =  cc*nvec(2) - a8*wpvec(2)*nvec(3)
  evr(3,3) = -cc*nvec(1) - a8*wpvec(3)*nvec(3)
  evr(4,3) = -a8*wpvec(4)*nvec(3)
  evr(5,3) =  cc*(wpvec(2)*nvec(2)-wpvec(3)*nvec(1)) - a10*nvec(3)
  evr(1,4) = h3
  evr(2,4) = wpvec(2) + h1*nvec(1)
  evr(3,4) = wpvec(3) + h1*nvec(2)
  evr(4,4) = wpvec(4) + h1*nvec(3)
  evr(5,4) = H + h1*V
  evr(1,5) = h3
  evr(2,5) = wpvec(2) + h2*nvec(1)
  evr(3,5) = wpvec(3) + h2*nvec(2)
  evr(4,5) = wpvec(4) + h2*nvec(3)
  evr(5,5) = H + h2*V
  evr(1,1) = evr(1,1)*h3
  evr(1,2) = evr(1,2)*h3
  evr(1,3) = evr(1,3)*h3
  evr(2,1) = evr(2,1)*h3
  evr(2,2) = evr(2,2)*h3
  evr(2,3) = evr(2,3)*h3
  evr(2,4) = evr(2,4)*h3
  evr(2,5) = evr(2,5)*h3
  evr(3,1) = evr(3,1)*h3
  evr(3,2) = evr(3,2)*h3
  evr(3,3) = evr(3,3)*h3
  evr(3,4) = evr(3,4)*h3
  evr(3,5) = evr(3,5)*h3
  evr(4,1) = evr(4,1)*h3
  evr(4,2) = evr(4,2)*h3
  evr(4,3) = evr(4,3)*h3
  evr(4,4) = evr(4,4)*h3
  evr(4,5) = evr(4,5)*h3
  evr(5,1) = evr(5,1)*h3
  evr(5,2) = evr(5,2)*h3
  evr(5,3) = evr(5,3)*h3
  evr(5,4) = evr(5,4)*h3
  evr(5,5) = evr(5,5)*h3
end subroutine RightEigenvec
! =============================================================================
!> Computes matrix product (T_p^-1) * (P^-1).
!!
!! @param wvec   vector of conservative variables
!! @param wpvec  vector of primitive variables
!! @param nvec   components of the unit normal vector
!! @param V      contravariant velocity
!! @param H      total enthalpy
!! @param theta  preconditioning parameter
!! @param rhop   derivative of density wrp. to pressure
!! @param rhoT   derivative of density wrp. to temperature
!! @param hp     derivative of enthalpy wrp. to pressure
!! @param hT     derivative of enthalpy wrp. to temperature
!! @param q2     total velocity squared
!! @param mat    resulting matrix
!!
subroutine MatprodTp1_P1( wvec,wpvec,nvec,V,H,theta,rhop,rhoT,hp,hT,q2,mat)
  use ModDataTypes
  implicit none
! parameters
  real(rtype), intent( in) :: V, H, theta, rhop, rhoT, hp, hT, q2
  real(rtype), intent( in) :: wvec(5), wpvec(5), nvec(3)
  real(rtype), intent(out) :: mat(5,5)
! local variables
  real(rtype) :: a1, a1g, ra1g, a4, a5, a5rp, a5c, a5c5, a5rt2, a7, cc, &
                 a14, a15, a15rt, a16, a16rt, a17, a17rt, h0, rhoT2, vc
! *****************************************************************************
  a1    = wvec(1)*rhop*hT + rhoT*(1.D0-wvec(1)*hp)
  a1g   = wvec(1)*theta*hT + rhoT*(1.D0-wvec(1)*hp)
  ra1g  = 1.D0/a1g
  a4    = a1*ra1g
  a5    = wvec(1)*hT*ra1g
  cc    = 0.5D0*Sqrt((V*V)*(a4-1.D0)*(a4-1.D0)+4.D0*a5)
  a7    = V*(a4-1.D0)/(4.D0*cc)
  h0    = (a5*wvec(1))/(a1*wpvec(5))
  a5rp  = rhop*h0
  a5c   = a5/cc
  a5c5  = 0.5D0*a5c
  rhoT2 = (rhoT*rhoT)/(a1*a1g)
  a5rt2 = a5rp - rhoT2
  a14   = h0*(1.D0-rhop*(H-q2)-wvec(1)*hp)
  a15   = (H-q2)*rhoT + wvec(1)*hT
  a15rt = rhoT*a15/(a1*a1g)
  a16   = (0.5D0+a7)/a1
  a16rt = a16*rhoT
  a17   = (0.5D0-a7)/a1
  a17rt = a17*rhoT
  vc    = nvec(1)*wpvec(2) + nvec(2)*wpvec(3) + nvec(3)*wpvec(4)
  mat(1,1) = a5c*(nvec(2)*wpvec(4)-nvec(3)*wpvec(3)) - nvec(1)*(a14+a15rt)
  mat(1,2) = (wpvec(2)*nvec(1))*a5rt2
  mat(1,3) = (wpvec(3)*nvec(1))*a5rt2 + a5c*nvec(3)
  mat(1,4) = (wpvec(4)*nvec(1))*a5rt2 - a5c*nvec(2)
  mat(1,5) = -nvec(1)*a5rt2
  mat(2,1) = a5c*(nvec(3)*wpvec(2)-nvec(1)*wpvec(4)) - nvec(2)*(a14+a15rt)
  mat(2,2) = (wpvec(2)*nvec(2))*a5rt2 - a5c*nvec(3)
  mat(2,3) = (wpvec(3)*nvec(2))*a5rt2
  mat(2,4) = (wpvec(4)*nvec(2))*a5rt2 + a5c*nvec(1)
  mat(2,5) = -nvec(2)*a5rt2
  mat(3,1) = a5c*(nvec(1)*wpvec(3)-nvec(2)*wpvec(2)) - nvec(3)*(a14+a15rt)
  mat(3,2) = (wpvec(2)*nvec(3))*a5rt2 + a5c*nvec(2)
  mat(3,3) = (wpvec(3)*nvec(3))*a5rt2 - a5c*nvec(1)
  mat(3,4) = (wpvec(4)*nvec(3))*a5rt2
  mat(3,5) = -nvec(3)*a5rt2
  mat(4,1) = a15*a16 - a5c5*vc
  mat(4,2) = a16rt*wpvec(2) + a5c5*nvec(1)
  mat(4,3) = a16rt*wpvec(3) + a5c5*nvec(2)
  mat(4,4) = a16rt*wpvec(4) + a5c5*nvec(3)
  mat(4,5) = -a16rt
  mat(5,1) = a15*a17 + a5c5*vc
  mat(5,2) = a17rt*wpvec(2) - a5c5*nvec(1)
  mat(5,3) = a17rt*wpvec(3) - a5c5*nvec(2)
  mat(5,4) = a17rt*wpvec(4) - a5c5*nvec(3)
  mat(5,5) = -a17rt
end subroutine MatprodTp1_P1
! =============================================================================
!> Computes matrix times the inverse of a similar matrix, where both matrices
!! have the structure of the preconditioning (or the transformation) matrix.
!! Thus, products like P*G^-1 or G*P^-1 can be evaluated efficiently.
!!
!! @param wpvec  vector of primitive variables
!! @param q2     total velocity squared
!! @param amat   matrix A
!! @param bmat   matrix B
!! @param cmat   resulting matrix C, i.e., C=A*B
!!
subroutine MatrixTimesInverse( wpvec,q2,amat,bmat,cmat)
  use ModDataTypes
  implicit none
! parameters
  real(rtype), intent( in) :: q2
  real(rtype), intent( in) :: wpvec(5), amat(5,5), bmat(5,5)
  real(rtype), intent(out) :: cmat(5,5)
! local variables
  real(rtype) :: a1, a2, a3, a4, a5
! *****************************************************************************
  a1 = amat(1,1)*bmat(1,1) + amat(1,5)*bmat(5,1)
  a2 = amat(1,1)*bmat(1,2) + amat(1,5)*bmat(5,2)
  a3 = amat(1,1)*bmat(1,3) + amat(1,5)*bmat(5,3)
  a4 = amat(1,1)*bmat(1,4) + amat(1,5)*bmat(5,4)
  a5 = amat(1,1)*bmat(1,5) + amat(1,5)*bmat(5,5)
  cmat(1,1) = a1
  cmat(2,1) = wpvec(2)*a1 - wpvec(2)
  cmat(3,1) = wpvec(3)*a1 - wpvec(3)
  cmat(4,1) = wpvec(4)*a1 - wpvec(4)
  cmat(5,1) = amat(5,1)*bmat(1,1) + amat(5,5)*bmat(5,1) - q2
  cmat(1,2) = a2
  cmat(2,2) = wpvec(2)*a2 + 1.D0
  cmat(3,2) = wpvec(3)*a2
  cmat(4,2) = wpvec(4)*a2
  cmat(5,2) = amat(5,1)*bmat(1,2) + amat(5,5)*bmat(5,2) + wpvec(2)
  cmat(1,3) = a3
  cmat(2,3) = wpvec(2)*a3
  cmat(3,3) = wpvec(3)*a3 + 1.D0
  cmat(4,3) = wpvec(4)*a3
  cmat(5,3) = amat(5,1)*bmat(1,3) + amat(5,5)*bmat(5,3) + wpvec(3)
  cmat(1,4) = a4
  cmat(2,4) = wpvec(2)*a4
  cmat(3,4) = wpvec(3)*a4
  cmat(4,4) = wpvec(4)*a4 + 1.D0
  cmat(5,4) = amat(5,1)*bmat(1,4) + amat(5,5)*bmat(5,4) + wpvec(4)
  cmat(1,5) = a5
  cmat(2,5) = wpvec(2)*a5
  cmat(3,5) = wpvec(3)*a5
  cmat(4,5) = wpvec(4)*a5
  cmat(5,5) = amat(5,1)*bmat(1,5) + amat(5,5)*bmat(5,5)
end subroutine MatrixTimesInverse
! =============================================================================
!> Computes matrix times vector (n=5).
!!
!! @param a  5x5 matrix A
!! @param v  vector v (length 5)
!! @param c  resulting vector c, i.e., c=A*v
!!
subroutine MatVecProd5( a,v,c)
  use ModDataTypes
  implicit none
! parameters
  real(rtype), intent( in) :: a(5,5), v(5)
  real(rtype), intent(out) :: c(5)
! *****************************************************************************
  c(1) = a(1,1)*v(1) + a(1,2)*v(2) + a(1,3)*v(3) + a(1,4)*v(4) + a(1,5)*v(5)
  c(2) = a(2,1)*v(1) + a(2,2)*v(2) + a(2,3)*v(3) + a(2,4)*v(4) + a(2,5)*v(5)
  c(3) = a(3,1)*v(1) + a(3,2)*v(2) + a(3,3)*v(3) + a(3,4)*v(4) + a(3,5)*v(5)
  c(4) = a(4,1)*v(1) + a(4,2)*v(2) + a(4,3)*v(3) + a(4,4)*v(4) + a(4,5)*v(5)
  c(5) = a(5,1)*v(1) + a(5,2)*v(2) + a(5,3)*v(3) + a(5,4)*v(4) + a(5,5)*v(5)
end subroutine MatVecProd5
!> Applies inflow boundary condition at dummy points (subsonic flow assumed).
!!
!! @param ibegn  indirect pointer to first node of the boundary
!! @param iendn  indirect pointer to last node of the boundary
!!
subroutine BcondInflow( ibegn,iendn)
  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModInterfaces, only : DependentVarsOne
  implicit none
! parameters
  integer, intent(in) :: ibegn, iendn
! local variables
  integer     :: ib, ibn, idn, ie
  real(rtype) :: ds, sxn, syn, rrho, u, v, uabs, unorm, cosa, c02, rinv, &
                 dis, cb, cc02, tb, pb, rhob, uabsb, ub, vb, gam1, ggm1, &
                 rgas
! *****************************************************************************
  do ib=ibegn,iendn
    ibn = bnode(1,ib)         ! boundary node
    idn = bnode(2,ib)         ! dummy node
    ie  = bnode(3,ib)         ! edge to dummy node
    ds  = Sqrt(sij(1,ie)**2+sij(2,ie)**2)
    sxn = sij(1,ie)/ds
    syn = sij(2,ie)/ds
    gam1  = dv(4,ibn) - 1.D0
    ggm1  = dv(4,ibn)/gam1
    rgas  = gam1*dv(5,ibn)/dv(4,ibn)
    rrho  = 1.D0/cv(1,ibn)
    u     = cv(2,ibn)*rrho
    v     = cv(3,ibn)*rrho
    uabs  = Sqrt(u*u+v*v)
    unorm = u*sxn + v*syn
    if(uabs < 1.D-20)then
      cosa = 1.D0
    else
      cosa = -unorm/uabs
    endif
    c02  = dv(3,ibn)*dv(3,ibn) + 0.5D0*gam1*uabs*uabs
    rinv = unorm - 2.D0*dv(3,ibn)/gam1
    dis  = (gam1*cosa*cosa+2.D0)*c02/(gam1*rinv*rinv) - 0.5D0*gam1
    if(dis < 0D0)then
      write(*,1000) ibn,dis
      dis = 1.D-20
    endif
    cb    = -rinv*(gam1/(gam1*cosa*cosa+2.D0))*(1.D0+cosa*Sqrt(dis))
    cc02  = Min(cb*cb/c02, 1.D0)
    tb    = ttinl*cc02
    pb    = ptinl*(tb/ttinl)**ggm1
    rhob  = pb/(rgas*tb)
    uabsb = 2.D0*dv(5,ibn)*(ttinl-tb)
    uabsb = Sqrt(uabsb)
    ub    = uabsb*Cos(betainl)
    vb    = uabsb*Sin(betainl)
    cv(1,idn) = rhob
    cv(2,idn) = rhob*ub
    cv(3,idn) = rhob*vb
    cv(4,idn) = pb/gam1 + 0.5D0*rhob*(ub*ub+vb*vb)
    call DependentVarsOne(idn)
  enddo  ! ib
1000  format(" Warning (BcondInflow): discriminant<0 at boundary node ",I3, &
             ", d= ",1PE13.5)
end subroutine BcondInflow
!> Prints user-input parameters for checking purposes.
!!
subroutine PrintParams
  use ModControl
  use ModFiles
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  implicit none
! local variables
  integer :: i
! *****************************************************************************
  write(*,1000) Trim( title)
! physics - general
  write(*,1005)
  write(*,1040) kflow,"E=external flow, I=internal flow"
  write(*,1040) kequs,"E=Euler, N=Navier-Stokes"
  write(*,1045) gamma,"ratio of specific heats"
  write(*,1045) cpgas,"specific heat coefficient (p=const.)"
  write(*,1045) renum,"Reynolds number"
  write(*,1045) refvel,"reference velocity"
  write(*,1045) refrho,"reference density"
  write(*,1045) refvisc,"laminar viscosity"
  write(*,1045) prlam,"laminar Prandtl number"
! physics - external/internal flow
  if(kflow=="e" .or. kflow=="E")then
    write(*,1010)
    write(*,1045) machinf,"Mach-number at infinity"
    write(*,1045) alpha*rad,"angle of attack [deg]"
    write(*,1045) pinf,"static pressure at infinity [Pa]"
    write(*,1045) tinf,"static temperature at infinity [K]"
  else
    write(*,1015)
    write(*,1045) ptinl,"total pressure at inlet [Pa]"
    write(*,1045) ttinl,"total temperature at inlet [K]"
    write(*,1045) betainl*rad,"flow angle at inlet (with x-axis) [deg]"
    write(*,1045) pout,"static pressure at outlet [Pa]"
    write(*,1045) betaout*rad,"approx. flow angle at outlet (with x-axis) [deg]"
    write(*,1045) p12rat,"approx. ratio of inlet to outlet static pressure"
  endif
! geometrical reference values
  write(*,1020)
  write(*,1045) xref,"x-coordinate of reference point (moment coefficient) [m]"
  write(*,1045) yref,"y-coordinate             - '' -"
  write(*,1045) cref,"reference or cord length [m]"
! iteration control
  write(*,1025)
  write(*,1050) maxiter,"max. number of iterations"
  write(*,1050) outstep,"number of iterations between solution dumps"
  write(*,1045) convtol,"2-norm of density change to stop the iteration"
  write(*,1040) lrest,"use previous solution for restart (Y=yes, N=no)"
! numerical parameters
  write(*,1030)
  write(*,1045) cfl,"CFL-number"
  write(*,1045) epsirs,"coefficient of implicit residual smoothing (<=0 - no smoothing)"
  write(*,1050) nitirs,"no. of Jacobi iterations for residual smoothing"
  write(*,1040) ktimst,"L=local, G=global time-stepping"
  write(*,1040) kprecond,"low Mach-number preconditioning (Y/N)"
  write(*,1045) precoeff,"preconditioning parameter K"
  write(*,1055) iorder,"1st-order (1) / 2nd-order (2) Roe scheme"
  write(*,1045) limfac,"limiter coefficient (Roe scheme)"
  write(*,1045) epsentr,"entropy correction coefficient (Roe scheme)"
  write(*,1040) lvort,"correction of far-field due to single vortex (external flow)"
  write(*,1055) nrk,"number of Runge-Kutta stages (steady flow only)"
  write(*,1060) (ark  (i), i=1,nrk)
  write(*,1060) (betrk(i), i=1,nrk)
  write(*,1065) (ldiss(i), i=1,nrk)
! quantities to plot
  write(*,1035)
  do i=1,mxquant
    write(*,1041) lquant(i),Trim( cquant(i))
  enddo
  write(*,1070)
! formats
1000  format(A)
1005  format("#",/,"# Physics - general",/,"# ",17("-"))
1010  format("#",/,"# Physics - external flow",/,"# ",23("-"))
1015  format("#",/,"# Physics - internal flow",/,"# ",23("-"))
1020  format("#",/,"# Geometrical reference values",/,"# ",28("-"))
1025  format("#",/,"# Iteration control",/,"# ",17("-"))
1030  format("#",/,"# Numerical parameters",/,"# ",20("-"))
1035  format("#",/,"# Quantities to plot",/,"# ",18("-"))
1040  format(2X,A1,12X,"# ",A)
1041  format(2X,A1,2X,"# ",A)
1045  format(1X,1PE11.4,3X,"# ",A)
1050  format(2X,I6,7X,"# ",A)
1055  format(2X,I1,12X,"# ",A)
1060  format(5(2X,F6.4))
1065  format(5(2X,I6))
1070  format(80("-"))
end subroutine PrintParams
!> Applies far-field boundary condition to dummy points. Characteristic
!! boundary conditions are employed in the case of subsonic flow. Conservative
!! variables are inter- / extrapolated in the case of supersonic flow. Vortex
!! correction is optionally applied to the flow variables (subsonic only).
!!
!! @param ibegn  indirect pointer to first node of the boundary
!! @param iendn  indirect pointer to last node of the boundary
!! @param rhof   work array for density at the boundary
!! @param uf     work array for u-velocity at the boundary
!! @param vf     work array for v-velocity at the boundary
!! @param pf     work array for pressure at the boundary
!!
subroutine BcondFarfield( ibegn,iendn,rhof,uf,vf,pf)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : DependentVarsOne, CompTheta, Forces
  implicit none
! parameters
  integer, intent(in) :: ibegn, iendn
  real(rtype) :: rhof(:), uf(:), vf(:), pf(:)
! local variables
  integer     :: ib, ibn, idn, ie
  real(rtype) :: gmr, gmg, bet, cir, xa, ya, dist, angle, sn, dn, vc, qv2
  real(rtype) :: ds, sxn, syn, rhoe, ue, ve, qqe, pe, qn, crho0, &
                 rhoa, ua, va, pa, ul, vl, pl, sgn, pb, gam1, ggm1
  real(rtype) :: rhop, rhoT, hT, theta, a1, ra1g, a4, a5, cs
! *****************************************************************************
! free-stream values (optionally corrected by a vortex) -----------------------
! values corrected
  if(lvort == "Y")then
    call Forces
    bet = Sqrt(1.D0-machinf*machinf)
    cir = 0.25D0*cref*cl*qinf/pi
    do ib=ibegn,iendn
      ibn      = bnode(1,ib)         ! boundary node
      gam1     = dv(4,ibn) - 1.D0
      ggm1     = dv(4,ibn)/gam1
      gmr      = 1.D0/dv(4,ibn)
      gmg      = gam1/dv(4,ibn)
      xa       = x(ibn) - xref
      ya       = y(ibn) - yref
      dist     = Sqrt(xa*xa+ya*ya)
      angle    = Atan2(ya,xa)
      sn       = Sin(angle-alpha)
      dn       = 1.D0 - machinf*machinf*sn*sn
      vc       = cir*bet/(dn*dist)
      uf(ib)   = uinf + vc*Sin(angle)
      vf(ib)   = vinf - vc*Cos(angle)
      qv2      = uf(ib)*uf(ib) + vf(ib)*vf(ib)
      pf(ib)   = (pinf**gmg+gmg*rhoinf*(qinf*qinf-qv2)/(2.D0*pinf**gmr))**ggm1
      rhof(ib) = rhoinf*(pf(ib)/pinf)**gmr
    enddo
! not corrected
  else
    do ib=ibegn,iendn
      rhof(ib) = rhoinf
      uf(ib)   = uinf
      vf(ib)   = vinf
      pf(ib)   = pinf
    enddo
  endif
! computation of the boundary values ------------------------------------------
  do ib=ibegn,iendn
    ibn = bnode(1,ib)         ! boundary node
    idn = bnode(2,ib)         ! dummy node
    ie  = bnode(3,ib)         ! edge to dummy node
    ds  = Sqrt(sij(1,ie)**2+sij(2,ie)**2)
    sxn = sij(1,ie)/ds
    syn = sij(2,ie)/ds
    gam1  = dv(4,ibn) - 1.D0
    rhoe  = cv(1,ibn)
    ue    = cv(2,ibn)/rhoe
    ve    = cv(3,ibn)/rhoe
    qqe   = ue*ue + ve*ve
    pe    = dv(1,ibn)
    if(machinf < 1.D0)then
! --- subsonic flow (qn<0: inflow / qn>0: outflow)
      qn = sxn*ue + syn*ve
      if(kprecond == "Y")then
        rhop  =  rhoe/pe
        rhoT  = -rhoe/dv(2,ibn)
        hT    =  dv(5,ibn)
        theta = CompTheta( dv(4,ibn),dv(3,ibn),qqe)
        a1    = rhoe*rhop*hT + rhoT
        ra1g  = 1.D0/(rhoe*theta*hT + rhoT)
        a4    = a1*ra1g
        a5    = rhoe*hT*ra1g
        cs    = Sqrt((qn*qn)*((a4-1.D0)*(a4-1.D0))+4.D0*a5)
        crho0 = rhoe*cs
      else
        crho0 = dv(3,ibn)*rhoe
      endif
      if(qn < 0D0)then
        rhoa = rhof(ib)
        ua   = uf(ib)
        va   = vf(ib)
        pa   = pf(ib)
        ul   = ue
        vl   = ve
        pl   = pe
        sgn  = -1.D0
        pb   = 0.5D0*(pa+pl-crho0*(sxn*(ua-ul)+syn*(va-vl)))
      else
        rhoa = rhoe
        ua   = ue
        va   = ve
        pa   = pe
        ul   = uf(ib)
        vl   = vf(ib)
        pl   = pf(ib)
        sgn  = +1.D0
        pb   = pf(ib)
      endif
      cv(1,idn) = rhoa + (pb-pa)/(dv(3,ibn)**2)
      cv(2,idn) = cv(1,idn)*(ua+sgn*sxn*(pa-pb)/crho0)
      cv(3,idn) = cv(1,idn)*(va+sgn*syn*(pa-pb)/crho0)
      cv(4,idn) = pb/gam1 + 0.5D0*(cv(2,idn)**2+cv(3,idn)**2)/cv(1,idn)
    else
! --- supersonic flow (qn<0: inflow / qn>0: outflow)
      qn = sxn*ue + syn*ve
      if(qn < 0D0)then
        cv(1,idn) = rhoinf
        cv(2,idn) = rhoinf*uinf
        cv(3,idn) = rhoinf*vinf
        cv(4,idn) = pinf/gam1 + 0.5D0*rhoinf*qinf*qinf
      else
        cv(1,idn) = rhoe
        cv(2,idn) = rhoe*ue
        cv(3,idn) = rhoe*ve
        cv(4,idn) = pe/gam1 + 0.5D0*rhoe*qqe
      endif
    endif
    call DependentVarsOne(idn)
  enddo  ! ib
end subroutine BcondFarfield
!> Reads in a single-character user option from file (in ASCII format).
!!
!! @param iunit input-file unit
!! @return option as single letter
!!
function ReadChar( iunit)
  use ModDataTypes
  implicit none
! parameters
  integer, intent(in) :: iunit
! result
  character(1) :: ReadChar
! local variables
  character(chrlen) :: str
  integer :: i
! *****************************************************************************
  read(iunit,"(A)") str
  do i=1,Len_trim(str)
    ReadChar = str(i:i)
    if(ReadChar /= " ") exit  ! first non-empty char should be the option ...
  enddo
end function ReadChar
!> Computes convective fluxes in the normal direction at solid walls.
!!
subroutine FluxWalls
  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModNumerics
  implicit none
! local variables
  integer     :: i, j, ib, ibf, ibegf, iendf
  real(rtype) :: sx, sy, pl, pr
! *****************************************************************************
  ibegf = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    if(btype(ib)>=300 .and. btype(ib)<500)then
      do ibf=ibegf,iendf
        i        = bface(1,ibf)
        j        = bface(2,ibf)
        sx       = sbf(1,ibf)/12.D0
        sy       = sbf(2,ibf)/12.D0
        pl       = 5.D0*dv(1,i) +      dv(1,j)
        pr       =      dv(1,i) + 5.D0*dv(1,j)
        rhs(2,i) = rhs(2,i) + sx*pl
        rhs(3,i) = rhs(3,i) + sy*pl
        rhs(2,j) = rhs(2,j) + sx*pr
        rhs(3,j) = rhs(3,j) + sy*pr
      enddo
    endif
    ibegf = iendf + 1
  enddo
end subroutine FluxWalls
!> Reads in grid data and boundary segments.
!!
subroutine ReadGrid
  use ModFiles
  use ModGeometry
  use ModInterfaces, only : DummyNodes, ErrorMessage
  implicit none
! local variables
  integer :: errFlag, i, ib, ibn, ibf, ibegf, iendf, ibegn, iendn
! *****************************************************************************
  open(unit=ifGrid, file=fnGrid, status="old", action="read", iostat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot open grid file")
  read(ifGrid,"(1X)")
  read(ifGrid,"(1X)")
  read(ifGrid,"(1X)")
! numbers of physical nodes, triangles and boundary segments
  read(ifGrid,*) nndint,ntria,nsegs
! boundary type, no. of boundary faces & nodes, boundary name
  allocate( btype(nsegs),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for btype()")
  allocate( bname(nsegs),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for bname()")
  allocate( ibound(2,nsegs),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for ibound()")
  read(ifGrid,"(1X)")
  do ib=1,nsegs
    read(ifGrid,  * ) btype(ib),ibound(1,ib),ibound(2,ib)
    read(ifGrid,"(A)") bname(ib)
  enddo
  nbfaces = ibound(1,nsegs)
  nbnodes = ibound(2,nsegs)
! definition of boundary faces / periodic nodes
  allocate( bnode(3,nbnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for bnode()")
  allocate( bface(2,nbfaces),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for bface()")
  do ibn=1,nbnodes
    bnode(1,ibn) = -777
    bnode(2,ibn) = -777      ! set in DummyNodes
    bnode(3,ibn) = -777      ! set in EdgesFinalize
  enddo
  do ibf=1,nbfaces
    bface(1,ibf) = -777
    bface(2,ibf) = -777
  enddo
  read(ifGrid,"(1X)")
  ibegf = 1
  ibegn = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    iendn = ibound(2,ib)
    if(btype(ib)>=700 .and. btype(ib)<800)then   ! periodic nodes
      do ibn=ibegn,iendn
        read(ifGrid,*) bnode(1,ibn),bnode(2,ibn)
      enddo
    else                                           ! boundary faces
      do ibf=ibegf,iendf
        read(ifGrid,*) bface(1,ibf),bface(2,ibf)
      enddo
    endif
    ibegf = iendf + 1
    ibegn = iendn + 1
  enddo
! check boundary faces pointer
  do ibf=1,nbfaces
    if(bface(1,ibf)<0 .or. bface(2,ibf)<0)then
      call ErrorMessage("array bface() not completely defined")
    endif
  enddo
! generate dummy nodes
  call DummyNodes
! grid nodes
  allocate( x(nnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for x()")
  allocate( y(nnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for y()")
  read(ifGrid,"(1X)")
  do i=1,nndint
    read(ifGrid,*) x(i),y(i)
  enddo
! triangles
  allocate( tria(3,ntria),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for tria()")
  read(ifGrid,"(1X)")
  do i=1,ntria
    read(ifGrid,*) tria(1,i),tria(2,i),tria(3,i)
  enddo
  close(unit=ifGrid)
end subroutine ReadGrid
!> Computes pressure forces and moments acting on the body. The contributions
!! are summed up by looping over ALL walls. Subroutine also computes lift and
!! drag coefficients.
!!
subroutine Forces
  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModPlotQuant
  implicit none
! local variables
  integer     :: ib, ibf, ibegf, iendf, n1, n2
  real(rtype) :: sx, sy, pwall, cp, xa, ya, dcx, dcy, cx, cy
! *****************************************************************************
! initialize force coefficients
  cx = 0D0
  cy = 0D0
  cm = 0D0
! loop over boundaries searching for walls
  ibegf = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    if(btype(ib)>=300 .and. btype(ib)<500)then
      do ibf=ibegf,iendf
        n1    = bface(1,ibf)
        n2    = bface(2,ibf)
        sx    = sbf(1,ibf)
        sy    = sbf(2,ibf)
        pwall = 0.5D0*(dv(1,n1)+dv(1,n2))
        cp    = 2.D0*(pwall-pinf)/(rhoinf*qinf*qinf)
        xa    = (0.5D0*(x(n1)+x(n2))-xref)/cref
        ya    = (0.5D0*(y(n1)+y(n2))-yref)/cref
        dcy   = sy*cp
        dcx   = sx*cp
        cy    = cy + dcy
        cx    = cx + dcx
        cm    = cm + dcx*ya - dcy*xa
      enddo
    endif ! btype
    ibegf = iendf + 1
  enddo   ! ib
! final lift and drag coefficients (pressure forces only!)
  cl = cy*Cos(alpha) - cx*Sin(alpha)
  cd = cy*Sin(alpha) + cx*Cos(alpha)
end subroutine Forces
!> Reads in user-input parameters.
!!
!! @param fname  path and name of the user input file
!!
subroutine ReadParams( fname)
  use ModControl
  use ModFiles
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : ErrorMessage, ReadChar
  implicit none
! parameters
  character(*), intent(in) :: fname
! local variables
  character(1) :: ch
  integer :: errFlag, i
! *****************************************************************************
  open(unit=ifInp, file=fname, status="old", action="read", iostat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot open input file")
  read(ifInp,"(1X)")
  read(ifInp,"(A)") title
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(A)") fnGrid
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(A)") fnFlow
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(A)") fnSurf
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(A)") fnConv
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(A)") fnRsti
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(A)") fnRsto
! physics - general
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  ch = ReadChar( ifInp)
  if(ch=="e" .or. ch=="E")then
    kflow = "E"
  else
    kflow = "I"
  endif
  ch = ReadChar( ifInp)
  if(ch=="e" .or. ch=="E")then
    kequs = "E"
  else
    kequs = "N"
  endif
  read(ifInp,*) gamma
  read(ifInp,*) cpgas
  read(ifInp,*) renum
  read(ifInp,*) refvel
  read(ifInp,*) refrho
  read(ifInp,*) prlam
! physics - external flow
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,*) machinf
  read(ifInp,*) alpha
  read(ifInp,*) pinf
  read(ifInp,*) tinf
! physics - internal flow
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,*) ptinl
  read(ifInp,*) ttinl
  read(ifInp,*) betainl
  read(ifInp,*) pout
  read(ifInp,*) betaout
  read(ifInp,*) p12rat
! geometrical reference values
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,*) xref
  read(ifInp,*) yref
  read(ifInp,*) cref
! iteration control
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,*) maxiter
  read(ifInp,*) outstep
  read(ifInp,*) convtol
  ch = ReadChar( ifInp)
  if(ch=="y" .or. ch=="Y")then
    lrest = "Y"
  else
    lrest = "N"
  endif
! numerical parameters
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,*) cfl
  read(ifInp,*) epsirs
  read(ifInp,*) nitirs
  ch = ReadChar( ifInp)
  if(ch=="l" .or. ch=="L")then
    ktimst = "L"
  else
    ktimst = "G"
  endif
  ch = ReadChar( ifInp)
  if(ch=="y" .or. ch=="Y")then
    kprecond = "Y"
  else
    kprecond = "N"
  endif
  read(ifInp,*) precoeff
  read(ifInp,*) iorder
  read(ifInp,*) limfac
  read(ifInp,*) epsentr
  ch = ReadChar( ifInp)
  if(ch=="y" .or. ch=="Y")then
    lvort = "Y"
  else
    lvort = "N"
  endif
  read(ifInp,*) nrk
  read(ifInp,*) (ark  (i), i=1,nrk)
  read(ifInp,*) (betrk(i), i=1,nrk)
  read(ifInp,*) (ldiss(i), i=1,nrk)
! quantities to plot
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  read(ifInp,"(1X)")
  do i=1,mxquant
    ch = ReadChar( ifInp)
    if(ch=="y" .or. ch=="Y")then
      lquant(i) = "Y"
    else
      lquant(i) = "N"
    endif
  enddo
! close input file
  close(unit=ifInp)
end subroutine ReadParams
!> Computes gradients of the density, u, v, and of the pressure with respect
!! to the x- and y-coordinates. Gradients are evaluated at the grid nodes.
!!
subroutine Gradients
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : Periodic
  implicit none
! local variables
  integer     :: i, j, ib, ibf, ibn, ie, ibegf, iendf, ibegn, iendn
  logical     :: flag
  real(rtype) :: rav, uav, vav, pav, sx, sy
  real(rtype) :: fcx(4), fcy(4)
! *****************************************************************************
! initialize gradients to zero
  do i=1,nnodes
    gradx(1,i) = 0D0
    gradx(2,i) = 0D0
    gradx(3,i) = 0D0
    gradx(4,i) = 0D0
    grady(1,i) = 0D0
    grady(2,i) = 0D0
    grady(3,i) = 0D0
    grady(4,i) = 0D0
  enddo
! sum up contributions at edge endpoints --------------------------------------
  do ie=1,nedint
    i = edge(1,ie)
    j = edge(2,ie)
! - average of variables
    rav = 0.5D0*(cv(1,i)+cv(1,j))
    uav = 0.5D0*(cv(2,i)/cv(1,i)+cv(2,j)/cv(1,j))
    vav = 0.5D0*(cv(3,i)/cv(1,i)+cv(3,j)/cv(1,j))
    pav = 0.5D0*(dv(1,i)+dv(1,j))
! - gradients (divided later by the volume)
    fcx(1) = rav*sij(1,ie)
    fcx(2) = uav*sij(1,ie)
    fcx(3) = vav*sij(1,ie)
    fcx(4) = pav*sij(1,ie)
    fcy(1) = rav*sij(2,ie)
    fcy(2) = uav*sij(2,ie)
    fcy(3) = vav*sij(2,ie)
    fcy(4) = pav*sij(2,ie)
    gradx(1,i) = gradx(1,i) + fcx(1)
    gradx(2,i) = gradx(2,i) + fcx(2)
    gradx(3,i) = gradx(3,i) + fcx(3)
    gradx(4,i) = gradx(4,i) + fcx(4)
    gradx(1,j) = gradx(1,j) - fcx(1)
    gradx(2,j) = gradx(2,j) - fcx(2)
    gradx(3,j) = gradx(3,j) - fcx(3)
    gradx(4,j) = gradx(4,j) - fcx(4)
    grady(1,i) = grady(1,i) + fcy(1)
    grady(2,i) = grady(2,i) + fcy(2)
    grady(3,i) = grady(3,i) + fcy(3)
    grady(4,i) = grady(4,i) + fcy(4)
    grady(1,j) = grady(1,j) - fcy(1)
    grady(2,j) = grady(2,j) - fcy(2)
    grady(3,j) = grady(3,j) - fcy(3)
    grady(4,j) = grady(4,j) - fcy(4)
  enddo
! contributions from the boundaries -------------------------------------------
  ibegf = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    flag  = .true.
    if(btype(ib)>=500 .and. btype(ib)<600) flag = .false.
    if(btype(ib)>=700 .and. btype(ib)<800) flag = .false.
    if(flag)then     ! all except symmetry and periodic boundaries
      do ibf=ibegf,iendf
        i   = bface(1,ibf)
        j   = bface(2,ibf)
        sx  = sbf(1,ibf)/12.D0
        sy  = sbf(2,ibf)/12.D0
! ----- node i
        rav = 5.D0*cv(1,i)         + cv(1,j)
        uav = 5.D0*cv(2,i)/cv(1,i) + cv(2,j)/cv(1,j)
        vav = 5.D0*cv(3,i)/cv(1,i) + cv(3,j)/cv(1,j)
        pav = 5.D0*dv(1,i)         + dv(1,j)
        gradx(1,i) = gradx(1,i) + rav*sx
        gradx(2,i) = gradx(2,i) + uav*sx
        gradx(3,i) = gradx(3,i) + vav*sx
        gradx(4,i) = gradx(4,i) + pav*sx
        grady(1,i) = grady(1,i) + rav*sy
        grady(2,i) = grady(2,i) + uav*sy
        grady(3,i) = grady(3,i) + vav*sy
        grady(4,i) = grady(4,i) + pav*sy
! ----- node j
        rav = 5.D0*cv(1,j)         + cv(1,i)
        uav = 5.D0*cv(2,j)/cv(1,j) + cv(2,i)/cv(1,i)
        vav = 5.D0*cv(3,j)/cv(1,j) + cv(3,i)/cv(1,i)
        pav = 5.D0*dv(1,j)         + dv(1,i)
        gradx(1,j) = gradx(1,j) + rav*sx
        gradx(2,j) = gradx(2,j) + uav*sx
        gradx(3,j) = gradx(3,j) + vav*sx
        gradx(4,j) = gradx(4,j) + pav*sx
        grady(1,j) = grady(1,j) + rav*sy
        grady(2,j) = grady(2,j) + uav*sy
        grady(3,j) = grady(3,j) + vav*sy
        grady(4,j) = grady(4,j) + pav*sy
      enddo
    endif
    ibegf = iendf + 1
  enddo
! correct at symmetry boundaries
  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if(btype(ib)>=500 .and. btype(ib)<600)then
      if(btype(ib)-500 < 2)then    ! x=const. line
        do ibn=ibegn,iendn
          i          = bnode(1,ibn)
          gradx(1,i) = 0D0
          grady(2,i) = 0D0
          gradx(3,i) = 0D0
          gradx(4,i) = 0D0
        enddo
      else                           ! y=const. line
        do ibn=ibegn,iendn
          i          = bnode(1,ibn)
          grady(1,i) = 0D0
          grady(2,i) = 0D0
          gradx(3,i) = 0D0
          grady(4,i) = 0D0
        enddo
      endif
    endif
    ibegn = iendn + 1
  enddo
! sum up at periodic boundaries
  call Periodic(gradx)
  call Periodic(grady)
! divide by the control volume ------------------------------------------------
  do i=1,nndint
    gradx(1,i) = gradx(1,i)/vol(i)
    gradx(2,i) = gradx(2,i)/vol(i)
    gradx(3,i) = gradx(3,i)/vol(i)
    gradx(4,i) = gradx(4,i)/vol(i)
    grady(1,i) = grady(1,i)/vol(i)
    grady(2,i) = grady(2,i)/vol(i)
    grady(3,i) = grady(3,i)/vol(i)
    grady(4,i) = grady(4,i)/vol(i)
  enddo
end subroutine Gradients
!> Reads in previous solution in order to restart the simulation. It also
!! reads the initial residual and number of previous iterations.
!!
subroutine ReadSolution
  use ModControl
  use ModFiles
  use ModGeometry
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : ErrorMessage
  implicit none
! local variables
  integer :: errFlag, i, n, nnodesDum, nconvDum
! *****************************************************************************
  open(unit=ifRsti, file=fnRsti, status="old", action="read", &
       form="unformatted", iostat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot open solution file")
! dimensions (for checking purposes)
  read(ifRsti) nnodesDum,nconvDum
  if(nnodesDum /= nnodes)  &
    call ErrorMessage("no. of nodes differs from the grid file")
  if(nconvDum /= nconv)  &
    call ErrorMessage("different number of conservative variables")
! initial residual, iteration # and solution
  read(ifRsti) drho1,iter
  read(ifRsti) ((cv(n,i), i=1,nnodes), n=1,nconv)
  close(ifRsti)
end subroutine ReadSolution
!> Computes gradients of the density, u, v, pressure and of the temperature
!! with respect to the x- and y-coordinates. Gradients are evaluated at
!! the grid nodes.
!!
subroutine GradientsVisc
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : Periodic
  implicit none
! local variables
  integer     :: i, j, ib, ibf, ibn, ie, ibegf, iendf, ibegn, iendn
  logical     :: flag
  real(rtype) :: rav, uav, vav, pav, tav, sx, sy
  real(rtype) :: fcx(5), fcy(5)
! *****************************************************************************
! initialize gradients to zero
  do i=1,nnodes
    gradx(1,i) = 0D0
    gradx(2,i) = 0D0
    gradx(3,i) = 0D0
    gradx(4,i) = 0D0
    gradx(5,i) = 0D0
    grady(1,i) = 0D0
    grady(2,i) = 0D0
    grady(3,i) = 0D0
    grady(4,i) = 0D0
    grady(5,i) = 0D0
  enddo
! sum up contributions at edge endpoints --------------------------------------
  do ie=1,nedint
    i = edge(1,ie)
    j = edge(2,ie)
! - average of variables
    rav = 0.5D0*(cv(1,i)+cv(1,j))
    uav = 0.5D0*(cv(2,i)/cv(1,i)+cv(2,j)/cv(1,j))
    vav = 0.5D0*(cv(3,i)/cv(1,i)+cv(3,j)/cv(1,j))
    pav = 0.5D0*(dv(1,i)+dv(1,j))
    tav = 0.5D0*(dv(2,i)+dv(2,j))
! - gradients (divided later by the volume)
    fcx(1) = rav*sij(1,ie)
    fcx(2) = uav*sij(1,ie)
    fcx(3) = vav*sij(1,ie)
    fcx(4) = pav*sij(1,ie)
    fcx(5) = tav*sij(1,ie)
    fcy(1) = rav*sij(2,ie)
    fcy(2) = uav*sij(2,ie)
    fcy(3) = vav*sij(2,ie)
    fcy(4) = pav*sij(2,ie)
    fcy(5) = tav*sij(2,ie)
    gradx(1,i) = gradx(1,i) + fcx(1)
    gradx(2,i) = gradx(2,i) + fcx(2)
    gradx(3,i) = gradx(3,i) + fcx(3)
    gradx(4,i) = gradx(4,i) + fcx(4)
    gradx(5,i) = gradx(5,i) + fcx(5)
    gradx(1,j) = gradx(1,j) - fcx(1)
    gradx(2,j) = gradx(2,j) - fcx(2)
    gradx(3,j) = gradx(3,j) - fcx(3)
    gradx(4,j) = gradx(4,j) - fcx(4)
    gradx(5,j) = gradx(5,j) - fcx(5)
    grady(1,i) = grady(1,i) + fcy(1)
    grady(2,i) = grady(2,i) + fcy(2)
    grady(3,i) = grady(3,i) + fcy(3)
    grady(4,i) = grady(4,i) + fcy(4)
    grady(5,i) = grady(5,i) + fcy(5)
    grady(1,j) = grady(1,j) - fcy(1)
    grady(2,j) = grady(2,j) - fcy(2)
    grady(3,j) = grady(3,j) - fcy(3)
    grady(4,j) = grady(4,j) - fcy(4)
    grady(5,j) = grady(5,j) - fcy(5)
  enddo
! contributions from the boundaries -------------------------------------------
  ibegf = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    flag  = .true.
    if(btype(ib)>=500 .and. btype(ib)<600) flag = .false.
    if(btype(ib)>=700 .and. btype(ib)<800) flag = .false.
    if(flag)then     ! all except symmetry and periodic boundaries
      do ibf=ibegf,iendf
        i   = bface(1,ibf)
        j   = bface(2,ibf)
        sx  = sbf(1,ibf)/12.D0
        sy  = sbf(2,ibf)/12.D0
! ----- node i
        rav = 5.D0*cv(1,i)         + cv(1,j)
        uav = 5.D0*cv(2,i)/cv(1,i) + cv(2,j)/cv(1,j)
        vav = 5.D0*cv(3,i)/cv(1,i) + cv(3,j)/cv(1,j)
        pav = 5.D0*dv(1,i)         + dv(1,j)
        tav = 5.D0*dv(2,i)         + dv(2,j)
        gradx(1,i) = gradx(1,i) + rav*sx
        gradx(2,i) = gradx(2,i) + uav*sx
        gradx(3,i) = gradx(3,i) + vav*sx
        gradx(4,i) = gradx(4,i) + pav*sx
        gradx(5,i) = gradx(5,i) + tav*sx
        grady(1,i) = grady(1,i) + rav*sy
        grady(2,i) = grady(2,i) + uav*sy
        grady(3,i) = grady(3,i) + vav*sy
        grady(4,i) = grady(4,i) + pav*sy
        grady(5,i) = grady(5,i) + tav*sy
! ----- node j
        rav = 5.D0*cv(1,j)         + cv(1,i)
        uav = 5.D0*cv(2,j)/cv(1,j) + cv(2,i)/cv(1,i)
        vav = 5.D0*cv(3,j)/cv(1,j) + cv(3,i)/cv(1,i)
        pav = 5.D0*dv(1,j)         + dv(1,i)
        tav = 5.D0*dv(2,j)         + dv(2,i)
        gradx(1,j) = gradx(1,j) + rav*sx
        gradx(2,j) = gradx(2,j) + uav*sx
        gradx(3,j) = gradx(3,j) + vav*sx
        gradx(4,j) = gradx(4,j) + pav*sx
        gradx(5,j) = gradx(5,j) + tav*sx
        grady(1,j) = grady(1,j) + rav*sy
        grady(2,j) = grady(2,j) + uav*sy
        grady(3,j) = grady(3,j) + vav*sy
        grady(4,j) = grady(4,j) + pav*sy
        grady(5,j) = grady(5,j) + tav*sy
      enddo
    endif
    ibegf = iendf + 1
  enddo
! correct at symmetry boundaries
  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if(btype(ib)>=500 .and. btype(ib)<600)then
      if(btype(ib)-500 < 2)then    ! x=const. line
        do ibn=ibegn,iendn
          i          = bnode(1,ibn)
          gradx(1,i) = 0D0
          grady(2,i) = 0D0
          gradx(3,i) = 0D0
          gradx(4,i) = 0D0
          gradx(5,i) = 0D0
        enddo
      else                           ! y=const. line
        do ibn=ibegn,iendn
          i          = bnode(1,ibn)
          grady(1,i) = 0D0
          grady(2,i) = 0D0
          gradx(3,i) = 0D0
          grady(4,i) = 0D0
          grady(5,i) = 0D0
        enddo
      endif
    endif
    ibegn = iendn + 1
  enddo
! sum up at periodic boundaries
  call Periodic(gradx)
  call Periodic(grady)
! divide by the control volume ------------------------------------------------
  do i=1,nndint
    gradx(1,i) = gradx(1,i)/vol(i)
    gradx(2,i) = gradx(2,i)/vol(i)
    gradx(3,i) = gradx(3,i)/vol(i)
    gradx(4,i) = gradx(4,i)/vol(i)
    gradx(5,i) = gradx(5,i)/vol(i)
    grady(1,i) = grady(1,i)/vol(i)
    grady(2,i) = grady(2,i)/vol(i)
    grady(3,i) = grady(3,i)/vol(i)
    grady(4,i) = grady(4,i)/vol(i)
    grady(5,i) = grady(5,i)/vol(i)
  enddo
end subroutine GradientsVisc
!> Integrates the four basic equations (continuity, momentum and energy) by
!! the explicit, multi-stage (Runge-Kutta) time-stepping scheme.
!!
!! @param iwork  integer work space for temporary variables
!! @param work   real work space for temporary variables
!!
subroutine Solver( iwork,work)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : BoundaryConditions, CompTheta, Cons2Prim, &
                            DependentVarsAll, DissipRoe1, DissipRoe1Prec, &
                            DissipRoe2, DissipRoe2Prec, ErrorMessage, &
                            FluxRoe1, FluxRoe2, FluxViscous, Gradients, &
                            GradientsVisc, Irsmoo, Limiter, LimiterInit, &
                            MatrixTimesInverse, Periodic, Prim2Cons, &
                            TimeStep, ZeroResiduals
  implicit none
! parameters
  integer     :: iwork(:)
  real(rtype) :: work(:)
! local variables
  integer     :: i, irk, mp, mp2
  real(rtype) :: blend1, fac, adtv, H, q2, rhop, rhoT, hT, theta, u, v
  real(rtype) :: wvec(5), wpvec(5), pmat(5,5), gmat1(5,5), dmat(5,5), r(5)
  real(rtype), allocatable :: dum1(:,:), dum2(:,:)
! *****************************************************************************
! calculate dimensions for dummy arrays (LimiterInit, Limiter, FluxRoe,
! Irsmoo, BoundaryConditions); check them
  mp  = nconv*nnodes
  mp2 = 2*mp
  if(mp2 > Ubound(work,1))then
    call ErrorMessage("insufficient work space in Solver")
  endif
! store previous solution; set dissipation = 0
  cvold(:,:) = cv(:,:)
  diss(:,:)  = 0D0
! compute the time step
  call TimeStep
! loop over the Runge-Kutta stages ============================================
  do irk=1,nrk
! - initialize dissipation
    if(irk>1 .and. ldiss(irk)/=0)then
      blend1 = 1.D0 - betrk(irk)
      do i=1,nnodes
        diss(1,i) = blend1*diss(1,i)
        diss(2,i) = blend1*diss(2,i)
        diss(3,i) = blend1*diss(3,i)
        diss(4,i) = blend1*diss(4,i)
      enddo
    endif
! - viscous flux (Navier-Stokes eqs.)
    if(ldiss(irk)/=0 .and. kequs=="N")then
      call GradientsVisc
      call FluxViscous(betrk(irk))
    endif
! - Roe's flux-difference splitting scheme (upwind)
    ! limiter and upwind dissipation
    if(ldiss(irk) /= 0)then
      if(iorder < 2)then
        if(kprecond == "Y")then
          call DissipRoe1Prec(betrk(irk))
        else
          call DissipRoe1(betrk(irk))
        endif
      else
!!!     dum1 = Reshape( work(1:mp)      ,(/4, nnodes/))
!!!     dum2 = Reshape( work((mp+1):mp2),(/4, nnodes/))
        allocate(dum1(4,nnodes))
        allocate(dum2(4,nnodes))
        if(kequs == "E")then
        call Gradients
        endif
        call LimiterInit(dum1,dum2)
        call Limiter(dum1,dum2)
        deallocate(dum1)
        deallocate(dum2)
        if(kprecond == "Y")then
          call DissipRoe2Prec(betrk(irk))
        else
          call DissipRoe2(betrk(irk))
        endif
      endif
    endif
    ! convective flux; add upwind dissipation => residual
    if(iorder < 2)then
      call FluxRoe1
    else
      call FluxRoe2
    endif
! - preconditioning
    if(kprecond == "Y")then
      do i=1,nndint
        rhop  =  cv(1,i)/dv(1,i)
        rhoT  = -cv(1,i)/dv(2,i)
        hT    = dv(5,i)
        u     = cv(2,i)/cv(1,i)
        v     = cv(3,i)/cv(1,i)
        q2    = u*u + v*v
        H     = (cv(4,i)+dv(1,i))/cv(1,i)
        theta = CompTheta( dv(4,i),dv(3,i),q2)
        wvec(1)  = cv(1,i)
        wvec(2)  = cv(2,i)
        wvec(3)  = cv(3,i)
        wvec(4)  = 0D0
        wvec(5)  = cv(4,i)
        wpvec(1) = dv(1,i)
        wpvec(2) = u
        wpvec(3) = v
        wpvec(4) = 0D0
        wpvec(5) = dv(2,i)
        call Cons2Prim(wvec,wpvec,H,q2,theta,rhoT,0D0,hT,gmat1)
        call Prim2Cons(wvec,wpvec,H,rhop,rhoT,0D0,hT,pmat)
        call MatrixTimesInverse(wpvec,q2,pmat,gmat1,dmat)
        r(1)     = rhs(1,i)
        r(2)     = rhs(2,i)
        r(3)     = rhs(3,i)
        r(4)     = rhs(4,i)
        rhs(1,i) = dmat(1,1)*r(1) + dmat(1,2)*r(2) + &
                   dmat(1,3)*r(3) + dmat(1,5)*r(4)
        rhs(2,i) = dmat(2,1)*r(1) + dmat(2,2)*r(2) + &
                   dmat(2,3)*r(3) + dmat(2,5)*r(4)
        rhs(3,i) = dmat(3,1)*r(1) + dmat(3,2)*r(2) + &
                   dmat(3,3)*r(3) + dmat(3,5)*r(4)
        rhs(4,i) = dmat(5,1)*r(1) + dmat(5,2)*r(2) + &
                   dmat(5,3)*r(3) + dmat(5,5)*r(4)
      enddo
    endif
! - correct residuals at symmetry/no-slip boundaries
    call ZeroResiduals
! - combine residuals at periodic boundaries
    call Periodic(rhs)
! - residual * time step / volume
    fac = ark(irk)*cfl
    do i=1,nndint
      adtv     = fac*tstep(i)/vol(i)
      rhs(1,i) = adtv*rhs(1,i)
      rhs(2,i) = adtv*rhs(2,i)
      rhs(3,i) = adtv*rhs(3,i)
      rhs(4,i) = adtv*rhs(4,i)
    enddo
! - implicit residual smoothing
    if(epsirs > 0D0)then
!!!   dum1 = Reshape( work(1:mp)      ,(/4, nnodes/))
!!!   dum2 = Reshape( work((mp+1):mp2),(/4, nnodes/))
      allocate(dum1(4,nnodes))
      allocate(dum2(4,nnodes))
      call Irsmoo(iwork,dum1,dum2)
      deallocate(dum1)
      deallocate(dum2)
      call ZeroResiduals
    endif
! - update - new solution, new dependent variables
    do i=1,nndint
      cv(1,i) = cvold(1,i) - rhs(1,i)
      cv(2,i) = cvold(2,i) - rhs(2,i)
      cv(3,i) = cvold(3,i) - rhs(3,i)
      cv(4,i) = cvold(4,i) - rhs(4,i)
    enddo
    call DependentVarsAll
! - boundary conditions
    call BoundaryConditions(work)
  enddo ! irk
end subroutine Solver
!> Initializes constants used by the solver.
!!
subroutine InitConstants
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none
! local variables
  real(rtype) :: gam1, rgas
! *****************************************************************************
  pi  = 4.D0*Atan(1.D0)
  rad = 180D0/pi
  if(kflow == 'E')then
! - external flow; it is assumed that "gamma" and "cpgas" specified in
!   the input file are valid for the complete far-field boundary
    gam1 = gamma - 1.D0
    rgas = gam1*cpgas/gamma
    alpha   = alpha/rad
    rhoinf  = pinf/(rgas*tinf)
    qinf    = machinf * Sqrt(gam1*cpgas*tinf)
    uinf    = qinf*Cos(alpha)
    vinf    = qinf*Sin(alpha)
    refrho  = rhoinf
    refvel  = qinf
    if(kequs == 'N')then
      refvisc = rhoinf*qinf*cref/renum
    else
      refvisc = 0D0
    endif
  else
! - internal flow
    betainl = betainl/rad
    betaout = betaout/rad
    if(kequs == 'N')then
      refvisc = refrho*refvel*cref/renum
    else
      refvisc = 0D0
    endif
  endif
end subroutine InitConstants
!> Computes maximum stable time step.
!!
subroutine TimeStep
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : CompTheta
  implicit none
! local variables
  integer     :: i
  real(rtype) :: sx, sy, ds, rrho, u, v, vc, cs, rhop, rhoT, hT, q2, &
                 theta, ra1g, a1, a4, a5, fmue, f1, f2, fac, dtv, cfac, &
                 lambdac, lambdav, tsmin
! *****************************************************************************
  if(kequs == "E")then
! - preconditioned Euler equations
    if(kprecond == "Y")then
      do i=1,nndint
        sx       = sproj(1,i)
        sy       = sproj(2,i)
        ds       = sx + sy
        u        = Abs(cv(2,i)/cv(1,i))
        v        = Abs(cv(3,i)/cv(1,i))
        rhop     = cv(1,i)/dv(1,i)
        rhoT     = -cv(1,i)/dv(2,i)
        hT       = dv(5,i)
        q2       = u*u + v*v
        theta    = CompTheta( dv(4,i),dv(3,i),q2)
        a1       = cv(1,i)*rhop*hT + rhoT
        ra1g     = 1.D0/(cv(1,i)*theta*hT+rhoT)
        a4       = a1*ra1g
        a5       = cv(1,i)*hT*ra1g
        vc       = (sx*u+sy*v)/ds
        cs       = Sqrt((vc*vc)*((a4-1.D0)*(a4-1.D0))+4.D0*a5)
        tstep(i) = (2.D0*vol(i))/((vc*(a4+1.D0)+cs)*ds)
      enddo
! - Euler equations
    else
      do i=1,nndint
        sx       = sproj(1,i)
        sy       = sproj(2,i)
        u        = Abs(cv(2,i)/cv(1,i))
        v        = Abs(cv(3,i)/cv(1,i))
        vc       = sx*u + sy*v
        cs       = dv(3,i)*(sx+sy)
        tstep(i) = vol(i)/(vc+cs)
      enddo
    endif
  else  ! kequs == "N"
    if(iorder > 1)then
      cfac = 1.D0
    else
      cfac = 2.D0
    endif
! - preconditioned Navier-Stokes equations (laminar)
    if(kprecond == "Y")then
      do i=1,nndint
        sx       = sproj(1,i)
        sy       = sproj(2,i)
        ds       = sx + sy
        rrho     = 1.D0/cv(1,i)
        u        = Abs(cv(2,i)*rrho)
        v        = Abs(cv(3,i)*rrho)
        rhop     = cv(1,i)/dv(1,i)
        rhoT     = -cv(1,i)/dv(2,i)
        hT       = dv(5,i)
        q2       = u*u + v*v
        theta    = CompTheta( dv(4,i),dv(3,i),q2)
        a1       = cv(1,i)*rhop*hT + rhoT
        ra1g     = 1.D0/(cv(1,i)*theta*hT+rhoT)
        a4       = a1*ra1g
        a5       = cv(1,i)*hT*ra1g
        vc       = (sx*u+sy*v)/ds
        cs       = Sqrt((vc*vc)*((a4-1.D0)*(a4-1.D0))+4.D0*a5)
        fmue     = dv(6,i)/prlam
        f1       = (4.D0*rrho)/3.D0
        f2       = dv(4,i)*rrho
        fac      = Max(f1,f2)
        dtv      = (fac*fmue)/vol(i)
        lambdac  = 0.5D0*((a4+1.D0)*vc+cs)*ds
        lambdav  = dtv*(sx*sx+sy*sy)
        tstep(i) = vol(i)/(lambdac+cfac*lambdav)
      enddo
! - Navier-Stokes equations (laminar)
    else
      do i=1,nndint
        sx       = sproj(1,i)
        sy       = sproj(2,i)
        rrho     = 1.D0/cv(1,i)
        u        = Abs(cv(2,i)*rrho)
        v        = Abs(cv(3,i)*rrho)
        fmue     = dv(6,i)/prlam
        f1       = (4.D0*rrho)/3.D0
        f2       = dv(4,i)*rrho
        fac      = Max(f1,f2)
        dtv      = (fac*fmue)/vol(i)
        vc       = sx*u + sy*v
        cs       = dv(3,i)*(sx+sy)
        lambdac  = vc + cs
        lambdav  = dtv*(sx*sx+sy*sy)
        tstep(i) = vol(i)/(lambdac+cfac*lambdav)
      enddo
    endif  ! kprecond
  endif    ! kequs
! in case of global time-stepping - find min. time step in domain -------------
  if(ktimst == "G")then
    tsmin = 1.D+32
    do i=1,nndint
      tsmin = Min(tsmin,tstep(i))
    enddo
    do i=1,nndint
      tstep(i) = tsmin
    enddo
  endif
end subroutine TimeStep
!> Initializes grid metrics: computes face vectors and cell volumes.
!!
!! @param niedge  pointer from a node to iedge()
!! @param iedge   linked list of edge endpoints:
!!                @li (1,*) = point j of edge (i,j)
!!                @li (2,*) = next point j which is also connected to i;
!!                            if <0 - no further connections
!!                @li (3,*) = pointer to edge() - used to associate face
!!                            vector sij() with the correct edge
!!
subroutine InitMetrics( niedge,iedge)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModInterfaces, only : ErrorMessage
  implicit none
! parameters
  integer, intent(in) :: niedge(:), iedge(:,:)
! local variables
  integer     :: cedge, d, errFlag, i, j, ic, ie, n
  real(rtype) :: x1, y1, x2, y2, x3, y3, area, pvol, cx, cy, sx, sy, vprod
! *****************************************************************************
! allocate memory
  allocate( sij(2,nedges),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for sij()")
  allocate( vol(nnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for vol()")
! zero out the variables
  do i=1,nnodes
    vol(i) = 0D0
  enddo
  do ie=1,nedges
    sij(1,ie) = 0D0
    sij(2,ie) = 0D0
  enddo
! loop over triangles, compute volumes and face vectors
  do ic=1,ntria
! - compute triangle area
    x1 = x(tria(1,ic))
    y1 = y(tria(1,ic))
    x2 = x(tria(2,ic))
    y2 = y(tria(2,ic))
    x3 = x(tria(3,ic))
    y3 = y(tria(3,ic))
    area = 0.5D0*((x1-x2)*(y1+y2)+(x2-x3)*(y2+y3)+(x3-x1)*(y3+y1))
    pvol = Abs(area)/3.D0
! - distribute area to the corner nodes (1/3 to each)
    vol(tria(1,ic)) = vol(tria(1,ic)) + pvol
    vol(tria(2,ic)) = vol(tria(2,ic)) + pvol
    vol(tria(3,ic)) = vol(tria(3,ic)) + pvol
! - compute center of the triangle
    cx = (x1+x2+x3)/3.D0
    cy = (y1+y2+y3)/3.D0
! - loop over individual nodes
    do n=1,3
      i = tria(n,ic)
      if(n < 3)then
        j = tria(n+1,ic)
      else
        j = tria(1,ic)
      endif
      if(i > j)then   ! lower index first
        d = i
        i = j
        j = d
      endif
! --- compute part of face vector associated with edge ij;
!     orient it to point in direction from i to j
      sx =  cy - 0.5D0*(y(i)+y(j))
      sy = -cx + 0.5D0*(x(i)+x(j))
      vprod = sx*(x(j)-x(i)) + sy*(y(j)-y(i))
      if(vprod < 0D0)then
        sx = -sx
        sy = -sy
      endif
! --- search corresponding edge and add (sx,sy) to sij()
      cedge = niedge(i)
10    continue
        if(iedge(1,cedge) == j)then
          ie        = iedge(3,cedge)
          sij(1,ie) = sij(1,ie) + sx
          sij(2,ie) = sij(2,ie) + sy
          goto 20
        endif
        cedge = iedge(2,cedge)
        if(cedge < 0)then
          call ErrorMessage("could not find edge to a node")
        endif
      goto 10
20    continue
    enddo   ! loop over nodes of a triangle
  enddo     ! loop over triangles
end subroutine InitMetrics
!> Computes projections of the control volumes on the x- and y-axis.
!! Variable (sproj) is used later to compute the time step.
!!
subroutine VolumeProjections
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModInterfaces, only : ErrorMessage
  implicit none
! local variables
  integer     :: errFlag, ibegn, iendn, ibegf, iendf
  integer     :: i, j, ib, ibf, ibn, ie
  real(rtype) :: sx, sy
! *****************************************************************************
! allocate memory
  allocate( sproj(2,nnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for sproj()")
! zero out the variable
  do i=1,nnodes
    sproj(1,i) = 0D0
    sproj(2,i) = 0D0
  enddo
! sum up contributions from interior edges
  do ie=1,nedint
    i          = edge(1,ie)
    j          = edge(2,ie)
    sx         = 0.5D0*Abs(sij(1,ie))
    sy         = 0.5D0*Abs(sij(2,ie))
    sproj(1,i) = sproj(1,i) + sx
    sproj(2,i) = sproj(2,i) + sy
    sproj(1,j) = sproj(1,j) + sx
    sproj(2,j) = sproj(2,j) + sy
  enddo
! add contributions from boundaries (except periodic)
  ibegf = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    if(btype(ib)<700 .or. btype(ib)>=800)then
      do ibf=ibegf,iendf       ! loop over boundary faces
        i          = bface(1,ibf)
        j          = bface(2,ibf)
        sx         = 0.25D0*Abs(sbf(1,ibf))
        sy         = 0.25D0*Abs(sbf(2,ibf))
        sproj(1,i) = sproj(1,i) + sx
        sproj(2,i) = sproj(2,i) + sy
        sproj(1,j) = sproj(1,j) + sx
        sproj(2,j) = sproj(2,j) + sy
      enddo
    endif
    ibegf = iendf + 1
  enddo
! sum up at periodic boundaries
  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if(btype(ib)>=700 .and. btype(ib)<800)then
      do ibn=ibegn,iendn
        i          = bnode(1,ibn)
        j          = bnode(2,ibn)
        sproj(1,i) = sproj(1,i) + sproj(1,j)
        sproj(1,j) = sproj(1,i)
        sproj(2,i) = sproj(2,i) + sproj(2,j)
        sproj(2,j) = sproj(2,i)
      enddo
    endif
    ibegn = iendn + 1
  enddo
end subroutine VolumeProjections
!> Duplicates control volumes at periodic boundaries. Initializes face
!! vectors at boundaries (except at periodic boundaries). Computes also
!! face vectors for dummy nodes (inlet/outlet/far-field). Sets coordinates
!! of dummy nodes equal to those of boundary nodes.
!!
!! @param marker  temporary node marker
!! @param btria   temporary pointer from boundary face to triangle
!!
subroutine InitMetricsBound( marker,btria)
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModInterfaces, only : ErrorMessage
  implicit none
! parameters
  integer :: marker(:), btria(:,:)
! local variables
  integer     :: d, errFlag, i, j, n1, n2, ibegf, ibegn, iendf, iendn, itype
  integer     :: ib, ibf, ibn, ic, ie, n
  logical     :: flag
  real(rtype) :: cx, cy, xm, ym, vprod
! *****************************************************************************
! combine control volumes at periodic nodes -----------------------------------
  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if(btype(ib)>=700 .and. btype(ib)<800)then
      do ibn=ibegn,iendn
        i      = bnode(1,ibn)
        j      = bnode(2,ibn)
        vol(i) = vol(i) + vol(j)
        vol(j) = vol(i)
      enddo
    endif
    ibegn = iendn + 1
  enddo
! compute face vectors at the boundaries --------------------------------------
! (find for each boundary face the corresp. triangle - required in order
! to define the face vectors such as to point OUT of the domain)
! reset node marker (nodes touched by bface()); mark nodes
  do i=1,nndint
    marker(i) = -1
  enddo
  do ibf=1,nbfaces
    marker(bface(1,ibf)) = 1
    marker(bface(2,ibf)) = 1
  enddo
! loop over triangles - if a face = bface()then store index of
! triangle in btria(1,*) => pointer from boundary face to triangle
  do ibf=1,nbfaces
    btria(1,ibf) = -1       ! empty now
  enddo
  do n=1,3                  ! loop over nodes of the triangles
    do ic=1,ntria
      i = tria(n,ic)
      if(n < 3)then       ! i, j define face of a triangle
        j = tria(n+1,ic)
      else
        j = tria(1,ic)
      endif
      if(i > j)then       ! lower index first
        d = i
        i = j
        j = d
      endif
      if(marker(i)==1 .and. marker(j)==1)then   ! must be on boundary
        do ibf=1,nbfaces                  ! search thru all boundary faces
          if(btria(1,ibf) == -1)then    ! triangle unknown
            n1 = bface(1,ibf)             ! n1, n2 define boundary face
            n2 = bface(2,ibf)
            if(n1 > n2)then             ! lower index first
              d  = n1
              n1 = n2
              n2 = d
            endif
            if(i==n1 .and. j==n2) btria(1,ibf) = ic  ! triangle found
          endif
        enddo
      endif
    enddo  ! loop over triangles
  enddo    ! loop over triangle nodes
! check if pointers from boundary faces to triangles valid
  do ibf=1,nbfaces
    if(btria(1,ibf) < 0)then
      call ErrorMessage("invalid pointer from boundary face to triangle")
    endif
  enddo
! allocate and compute boundary face vector sbf
  allocate( sbf(2,nbfaces),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for sbf()")
  do ibf=1,nbfaces
    sbf(1,ibf) = y(bface(2,ibf)) - y(bface(1,ibf))
    sbf(2,ibf) = x(bface(1,ibf)) - x(bface(2,ibf))
! - change orientation of sbf (must be outward facing)
    ic    = btria(1,ibf)
    cx    = (x(tria(1,ic))+x(tria(2,ic))+x(tria(3,ic)))/3.D0
    cy    = (y(tria(1,ic))+y(tria(2,ic))+y(tria(3,ic)))/3.D0
    xm    = 0.5D0*(x(bface(1,ibf))+x(bface(2,ibf)))
    ym    = 0.5D0*(y(bface(1,ibf))+y(bface(2,ibf)))
    vprod = sbf(1,ibf)*(xm-cx) + sbf(2,ibf)*(ym-cy)
    if(vprod < 0D0)then
      sbf(1,ibf) = -sbf(1,ibf)
      sbf(2,ibf) = -sbf(2,ibf)
    endif
  enddo
! generate face vectors for dummy edges ---------------------------------------
  ibegf = 1
  ibegn = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    iendn = ibound(2,ib)
    itype = btype(ib)
    flag  = .false.
    if(itype>=100 .and. itype<200) flag = .true.
    if(itype>=200 .and. itype<300) flag = .true.
    if(itype>=600 .and. itype<700) flag = .true.
    if(flag)then                      ! inlet/outlet/far-field boundary
      do n=1,2                          ! loop over nodes of boundary faces
        do ibf=ibegf,iendf              ! loop over boundary faces
          n1 = bface(n,ibf)             ! which node
          do i=ibegn,iendn              ! search for corresp. boundary node
            if(bnode(1,i) == n1)then
              ie        = bnode(3,i)    ! edge to dummy node
              sij(1,ie) = sij(1,ie) + 0.5D0*sbf(1,ibf)
              sij(2,ie) = sij(2,ie) + 0.5D0*sbf(2,ibf)
            endif
          enddo       ! i
        enddo         ! ibf
      enddo           ! n
    endif             ! flag
    ibegf = iendf + 1
    ibegn = iendn + 1
  enddo
! set coordinates and volumes of dummy nodes ----------------------------------
  do ie=nedint+1,nedges
    i      = edge(1,ie)        ! boundary node
    j      = edge(2,ie)        ! dummy node
    x(j)   = x(i)
    y(j)   = y(i)
    vol(j) = vol(i)
  enddo
end subroutine InitMetricsBound
!> Stores the current flow solution together with the initial residual
!! and the actual number of iterations.
!!
subroutine WriteSolution
  use ModControl
  use ModFiles
  use ModGeometry
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : ErrorMessage
  implicit none
! local variables
  integer :: errFlag, i, n
! *****************************************************************************
!!open(unit=ifRsto, file=fnRsto, status="unknown", action="write", &
!!     form="unformatted", iostat=errFlag)
!!if (errFlag /= 0) call ErrorMessage("cannot open solution file")
! dimensions
!!write(ifRsto) nnodes,nconv
! initial residual, iteration # and solution
!!write(ifRsto) drho1,iter
!!write(ifRsto) ((cv(n,i), i=1,nnodes), n=1,nconv)
!!close(ifRsto)
end subroutine WriteSolution
!> Allocates memory for conservative, dependent and for numerical variables.
!!
subroutine AllocateMemory
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : ErrorMessage
  implicit none
! local variables
  integer :: errFlag
! *****************************************************************************
! base flow variables
  allocate( cv(nconv,nnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for cv()")
  allocate( dv(ndepv,nnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for dv()")
! general numerical variables
  allocate( cvold(nconv,nnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for cvold()")
  allocate( diss(nconv,nnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for diss()")
  allocate( rhs(nconv,nnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for rhs()")
  allocate( tstep(nnodes),stat=errFlag)
  if(errFlag /= 0) call ErrorMessage("cannot allocate memory for tstep()")
! numerical variables required for higher-order Roe scheme
! and for viscous flow
  if(iorder > 1)then
    allocate( lim(nconv,nnodes),stat=errFlag)
    if(errFlag /= 0) call ErrorMessage("cannot allocate memory for lim()")
    if(kequs=="n" .or. kequs=="N")then
      allocate( gradx(5,nnodes),stat=errFlag)
      if(errFlag /= 0) call ErrorMessage("cannot allocate memory for gradx()")
      allocate( grady(5,nnodes),stat=errFlag)
      if(errFlag /= 0) call ErrorMessage("cannot allocate memory for grady()")
    else
      allocate( gradx(4,nnodes),stat=errFlag)
      if(errFlag /= 0) call ErrorMessage("cannot allocate memory for gradx()")
      allocate( grady(4,nnodes),stat=errFlag)
      if(errFlag /= 0) call ErrorMessage("cannot allocate memory for grady()")
    endif
  else
    if(kequs=="n" .or. kequs=="N")then
      allocate( gradx(5,nnodes),stat=errFlag)
      if(errFlag /= 0) call ErrorMessage("cannot allocate memory for gradx()")
      allocate( grady(5,nnodes),stat=errFlag)
      if(errFlag /= 0) call ErrorMessage("cannot allocate memory for grady()")
    endif
  endif
end subroutine AllocateMemory
!> Initializes conservative variables (initial guess). Ideal gas
!! is assumed (constant gas properties used everywhere).
!!
subroutine InitSolution
  use ModDataTypes
  use ModGeometry
  use ModPhysics
  implicit none
! local variables
  integer     :: i, j, ib, ibn, ibegn, iendn
  real(rtype) :: gam1, rgas, xmin, xmax, dx, pinl, temp, rho, &
                 cs, mach, q, dp, dbeta, beta, p, u, v
! *****************************************************************************
  gam1 = gamma - 1.D0
  rgas = gam1*cpgas/gamma
  if(kflow == "E")then
! - external flow
    do i=1,nnodes
      cv(1,i) = rhoinf
      cv(2,i) = rhoinf*uinf
      cv(3,i) = rhoinf*vinf
      cv(4,i) = pinf/gam1 + 0.5D0*rhoinf*qinf*qinf
    enddo
  else
! - internal flow; inlet assumed at xmin, outlet at xmax;
!   flow angle and pressure are linearly interpolated between
!   inlet and outlet
    xmin =  1.D+32
    xmax = -1.D+32
    do i=1,nndint
      xmin = Min(xmin,x(i))
      xmax = Max(xmax,x(i))
    enddo
    dx = xmax - xmin
    pinl  = p12rat*pout
    if(pinl >= ptinl) pinl = 0.99999D0*ptinl  ! otherwise reversed flow at inlet
    dp    = pout - pinl
    dbeta = betaout - betainl
    temp  = ttinl*(pinl/ptinl)**(gam1/gamma)
    rho   = pinl/(rgas*temp)
    cs    = Sqrt(gamma*pinl/rho)
    mach  = Sqrt(2.D0*((ttinl/temp)-1.D0)/gam1)
    q     = mach*cs
    do i=1,nnodes
      beta    = betainl + dbeta*(x(i)-xmin)/dx
      p       = pinl + dp*(x(i)-xmin)/dx
      rho     = p/(rgas*temp)
      u       = q*Cos(beta)
      v       = q*Sin(beta)
      cv(1,i) = rho
      cv(2,i) = rho*u
      cv(3,i) = rho*v
      cv(4,i) = p/gam1 + 0.5D0*rho*q*q
    enddo
  endif
! equalize flow variables at periodic nodes
  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if(btype(ib)>=700 .and. btype(ib)<800)then
      do ibn=ibegn,iendn
        i       = bnode(1,ibn)
        j       = bnode(2,ibn)
        cv(1,i) = 0.5D0*(cv(1,i)+cv(1,j))
        cv(2,i) = 0.5D0*(cv(2,i)+cv(2,j))
        cv(3,i) = 0.5D0*(cv(3,i)+cv(3,j))
        cv(4,i) = 0.5D0*(cv(4,i)+cv(4,j))
        cv(1,j) = cv(1,i)
        cv(2,j) = cv(2,i)
        cv(3,j) = cv(3,i)
        cv(4,j) = cv(4,i)
      enddo
    endif
    ibegn = iendn + 1
  enddo
end subroutine InitSolution
