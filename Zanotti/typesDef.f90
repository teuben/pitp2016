MODULE typesDef
  IMPLICIT NONE 
  PUBLIC 
  !
  ! ================================== This part of the typesDef can be modified by the user.  ================================== 
  !
  INTEGER, PARAMETER :: N = 2                               ! Polynomial degree of our approximation in space and time 
  INTEGER, PARAMETER :: nDim = 1                           ! The number of space dimensions that we actually want to simulate 
  REAL, PARAMETER    :: CFL = 0.2                           ! The Courant-Friedrichs-Lewy number < 1 

  INTEGER, PARAMETER :: nVar = 5                            ! The number of variables of the PDE system 

  !
  ! ==================================           Do NOT change the stuff below                 ==================================
  !
  ! The following variables contain important information about the numerical method. Do NOT change.  
  !
  INTEGER, PARAMETER :: d = 3                               ! This is the maximum number of space dimensions we want to deal with in our heads. !! NEVER change this parameter, unless you are bold and want to solve the Boltzmann equation !! 
  REAL, PARAMETER    :: PNPMTable(0:9) = (/ 1.0, 0.33, 0.17, 0.1, 0.069, 0.045,  0.038, 0.03, 0.02, 0.015 /)   ! maximum CFL numbers for PNPM schemes according to von Neumann stability analysis and experience     
  INTEGER            :: nDOF(0:3)                           ! The number of degrees of freedom in space and time 
  REAL               :: xiGPN(N+1), wGPN(N+1)               ! The Gauss-Legendre quadrature nodes and weights 
  REAL               :: xin(N+1)                            ! The nodes used for our basis (can in principle be different from the Gauss-Legendre nodes, like for example also the bad Gauss-Lobatto nodes) 
  REAL               :: MM(N+1,N+1), iMM(N+1,N+1)           ! Element mass-matrix and its inverse 
  REAL               :: Kxi(N+1,N+1), dudx(N+1,N+1)         ! Element stiffness matrix and discrete derivative operator 
  REAL               :: FLm(N+1,N+1), FLp(N+1,N+1)          ! Left flux matrices 
  REAL               :: FRm(N+1,N+1), FRp(N+1,N+1)          ! Right flux matrices 
  REAL               :: FLcoeff(N+1), FRcoeff(N+1)          ! extrapolation coefficients to the left and right boundary 
  REAL               :: F0(N+1), F1(N+1,N+1)                ! Time flux matrices 
  REAL               :: K1(N+1,N+1), iK1(N+1,N+1)           ! F1 - Ktau 
  INTEGER            :: dn(d)                               ! number of direct neighbors in each dimension 
  ! Stuff related to the problem setup, mesh and output 
  INTEGER            :: IMAX, JMAX, KMAX, NMAX              ! The number of cells in each space dimension & max. number of time steps 
  INTEGER            :: nElem, nFace, nNode                 ! The number of elements, faces and nodes 
  INTEGER            :: timestep                            ! the number of the current time step 
  REAL               :: xL(d), xR(d)                        ! computational domain 
  REAL               :: dx(d), dt                           ! The vector of mesh spacings in each dimension and the time step   
  REAL               :: time, tend                          ! current time  and final time 
  REAL               :: tio, dtio                           ! output time and output time interval 
  REAL, POINTER      :: x(:,:)                              ! the node coordinates (nDim, nNode) 
  INTEGER, PARAMETER :: nVtx = 2**nDim, nFac = 2*nDim       ! number of vertices and faces per element 
  INTEGER, POINTER   :: tri(:,:)                            ! connectivity from the element to the nodes 
  INTEGER, POINTER   :: Element2Face(:,:)                   ! connectivity from each element to its faces 
  CHARACTER(LEN=200) :: BaseFile                            ! Basic filename to write the results 
  REAL               :: SubOutputMatrix((N+1)**d,(N+1)**d)  ! Matrix needed for the plotting of the results on a fine subgrid 
  INTEGER            :: subtri(2**d,N**d)                   ! subcell connectivity (for fine output) 
  REAL               :: allsubxi(d,(N+1)**d)                ! subnodes (for fine output) 
  ! Some diagnostics                                        ! 
  REAL               :: tCPU1, tCPU2                        ! CPU times 
  INTEGER(8)         :: TEU                                 ! total element updates 
  !
  TYPE tFace
    REAL, POINTER    :: qL(:,:,:), qR(:,:,:)                ! pointer to left and right boundary-extrapolated state vector 
    REAL, POINTER    :: FL(:,:,:), FR(:,:,:)                ! pointer to left and right boundary-extrapolated flux vector 
    INTEGER          :: Left, Right                         ! pointer to left and right element 
    REAL             :: nv(d)                               ! face normal vector 
  END TYPE      
  TYPE(tFace), POINTER :: Face(:) 
  !
  ! The main variables of the ADER-DG scheme 
  !  
  REAL, POINTER      :: uh(:,:,:,:,:)                       ! the coefficients of the DG polynomial        (nVar, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: duh(:,:,:,:,:)                      ! the update coefficients of the DG polynomial (nVar, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: qhi(:,:,:,:,:)                      ! the time-averaged coefficients of the space-time predictor (nVar, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: Fhi(:,:,:,:,:,:)                    ! the time-averaged coefficients of the flux tensor of the space-time predictor (nVar, d, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: qbnd(:,:,:,:,:)                     ! the boundary-extrapolated data for the state vector Q in the element 
  REAL, POINTER      :: Fbnd(:,:,:,:,:)                     ! the boundary-extrapolated data for the normal flux F in the element 
  !
  ! Important info and parameters concerning the governing PDE system 
  !
  TYPE tEquations 
      REAL           :: gamma 
  END TYPE tEquations 
  
  TYPE(tEquations)   :: EQN 
  
  !
END MODULE typesDef 
    
    
    
    
