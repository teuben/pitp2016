SUBROUTINE ADERDGInit
    USE typesDef
    IMPLICIT NONE
    ! Local variables 
    INTEGER :: i, j, k, ii, jj, kk, l, c, iGP, iElem, VMAX(d), count, cnt  
    REAL    :: phi(N+1), phi0(N+1), phi1(N+1), phi_xi(N+1)
    REAL    :: phi_i(N+1), phi_j(N+1), phi_k(N+1) 
    REAL    :: u0(nVar), xGP(d), x0(d), subxi(N+1), aux(d), nv(d)
    INTEGER, POINTER :: idxn(:,:,:), idxe(:,:,:)  
    ! 
    ! ----------------- Some important preliminary stuff. Do not touch! ------------------- 
    dn(:) = 0 
    DO i = 1, nDim
        dn(i) = 1
    ENDDO 
    ! According to the number of space dimensions, we set the number of degrees of freedom 
    ! i = 0 is the time dimension 
    nDOF(:) = 1 
    DO i = 0, nDim
        nDOF(i) = N+1
    ENDDO 
    ! ------------------------------------------------------------------------------------- 
    !
    ! Some info about the PDE system 
    !
    EQN%gamma = 1.4                                                 ! ratio of specific heats for compressible Euler 
    !
    ! Here, you need to define the computational domain and the number of cells. 
    ! This is typically read from a parameter file 
    !
    xL = (/ 0.0 , -0.5, -0.5 /)                                     ! lower-left corner of the domain 
    xR = (/ 1.0,   0.5, +0.5 /)                                     ! upper right corner of the domain 
    IMAX = 50                                                       ! Number of elements in x,y,z direction 
    JMAX = 1 
    KMAX = 1  
    VMAX = (/ IMAX, JMAX, KMAX /)                                   ! Vector of the number of elements in each space dimension 
    dx = (xR-xL)/VMAX                                               ! Mesh spacing 
    NMAX = 100000                                                   ! Max. number of time steps 
    timestep = 0                                                    ! initial time step number 
    time = 0.                                                       ! initial time 
    tend = 0.2 !25                                                     ! final time 
    Basefile = 'Test'                                               ! Base filename for writing results 
    !
    nElem = IMAX*JMAX*KMAX                                          ! Number of elements 
    ALLOCATE(  uh(nVar, nDOF(1), nDOF(2), nDOF(3), nElem) )         ! Allocate the spatial degrees of freedom 
    ALLOCATE( duh(nVar, nDOF(1), nDOF(2), nDOF(3), nElem) )         ! Allocate the degrees of freedom for the update 
    ALLOCATE( qhi(nVar, nDOF(1), nDOF(2), nDOF(3), nElem) )         ! Allocate the time-averaged space-time degrees of freedom 
    ALLOCATE( Fhi(nVar, d, nDOF(1), nDOF(2), nDOF(3), nElem) )      ! Allocate the time-averaged space-time degrees of freedom 
    ALLOCATE( qbnd(nVar, 6, nDOF(2), nDOF(3), nElem) )              ! Allocate the time-averaged boundary-extrapolated values for Q 
    ALLOCATE( Fbnd(nVar, 6, nDOF(2), nDOF(3), nElem) )              ! Allocate the time-averaged boundary-extrapolated values for the normal flux F * n 
    nNode = (IMAX+1)*(JMAX+1)*(KMAX+1)                              ! number of nodes 
    ALLOCATE( x(d, nNode) )                                         ! Allocate the nodes 
    ALLOCATE( idxn(IMAX+dn(1),JMAX+dn(2),KMAX+dn(3))  )                                              
    ! Define the node coordinates and the node numbers                                 
    count = 0 
    DO k = 1, KMAX+dn(3)  
        DO j = 1, JMAX+dn(2) 
            DO i = 1, IMAX+dn(1) 
                count = count + 1 
                x(:,count) = xL(:) + (/ i-1, j-1, k-1/)*dx(:) 
                idxn(i,j,k) = count 
            ENDDO
        ENDDO
    ENDDO
    ! define the connectivity between the elements and the nodes. You can do this more compactly via a loop, but I prefer the select case, which makes the explicit 
    ! construction of each element clearer 
    ALLOCATE( tri(nVtx, nElem)     ) 
    ALLOCATE( idxe(IMAX,JMAX,KMAX) ) 
    count = 0  
    DO k = 1, KMAX  
        DO j = 1, JMAX
            DO i = 1, IMAX 
                count = count + 1
                idxe(i,j,k) = count 
                SELECT CASE(nDim)
                CASE(1)                    
                    tri(:,count) = (/ idxn(i,j,k), idxn(i+1,j,k) /) 
                CASE(2)
                    tri(:,count) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i,j+1,k), idxn(i+1,j+1,k) /) 
                CASE(3)
                    tri(:,count) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i,j+1,k), idxn(i+1,j+1,k), idxn(i,j,k+1), idxn(i+1,j,k+1), idxn(i,j+1,k+1), idxn(i+1,j+1,k+1) /) 
                END SELECT                
            ENDDO
        ENDDO
    ENDDO
    DEALLOCATE( idxn )
    ! define the connectivity between the faces and the elements 
    ! count how many faces we have 
    nFace = KMAX*JMAX*(IMAX+1)
    IF(nDim>=2) THEN 
        nFace = nFace + KMAX*(JMAX+1)*IMAX
    ENDIF
    IF(nDim>=3) THEN 
        nFace = nFace + (KMAX+1)*JMAX*IMAX 
    ENDIF
    ALLOCATE( Face(nFace) ) 
    ! x faces 
    count = 0 
    DO k = 1, KMAX
     DO j = 1, JMAX
      DO i = 1, IMAX+1  
          count = count + 1 
          IF(i.EQ.1) THEN
            Face(count)%Left  = 0 
            Face(count)%Right = idxe(i,j,k)
            ALLOCATE( Face(count)%qL(nVar,nDOF(2),nDOF(3)) ) 
            ALLOCATE( Face(count)%FL(nVar,nDOF(2),nDOF(3)) ) 
            Face(count)%qR => qBnd(:,1,:,:,Face(count)%Right) 
            Face(count)%FR => FBnd(:,1,:,:,Face(count)%Right) 
          ELSEIF(i.EQ.IMAX+1) THEN 
            Face(count)%Left  = idxe(i-1,j,k) 
            Face(count)%Right = 0 
            Face(count)%qL => qBnd(:,2,:,:,Face(count)%Left ) 
            Face(count)%FL => FBnd(:,2,:,:,Face(count)%Left ) 
            ALLOCATE( Face(count)%qR(nVar,nDOF(2),nDOF(3)) ) 
            ALLOCATE( Face(count)%FR(nVar,nDOF(2),nDOF(3)) ) 
          ELSE              
            Face(count)%Left  = idxe(i-1,j,k) 
            Face(count)%Right = idxe(i,j,k) 
            Face(count)%qL => qBnd(:,2,:,:,Face(count)%Left ) 
            Face(count)%qR => qBnd(:,1,:,:,Face(count)%Right) 
            Face(count)%FL => FBnd(:,2,:,:,Face(count)%Left ) 
            Face(count)%FR => FBnd(:,1,:,:,Face(count)%Right) 
          ENDIF 
          Face(count)%nv = (/ 1., 0., 0. /) ! set face normal vector 
      ENDDO
     ENDDO
    ENDDO 
    ! y faces 
    IF(nDim>=2) THEN
        DO k = 1, KMAX
         DO j = 1, JMAX+1 
          DO i = 1, IMAX  
              count = count + 1 
              IF(j.EQ.1) THEN
                Face(count)%Left  = 0 
                Face(count)%Right = idxe(i,j,k)
                ALLOCATE( Face(count)%qL(nVar,nDOF(2),nDOF(3)) ) 
                ALLOCATE( Face(count)%FL(nVar,nDOF(2),nDOF(3)) ) 
                Face(count)%qR => qBnd(:,3,:,:,Face(count)%Right) 
                Face(count)%FR => FBnd(:,3,:,:,Face(count)%Right)
              ELSEIF(j.EQ.JMAX+1) THEN 
                Face(count)%Left  = idxe(i,j-1,k) 
                Face(count)%Right = 0 
                Face(count)%qL => qBnd(:,4,:,:,Face(count)%Left ) 
                Face(count)%FL => FBnd(:,4,:,:,Face(count)%Left ) 
                ALLOCATE( Face(count)%qR(nVar,nDOF(2),nDOF(3)) ) 
                ALLOCATE( Face(count)%FR(nVar,nDOF(2),nDOF(3)) ) 
              ELSE              
                Face(count)%Left  = idxe(i,j-1,k) 
                Face(count)%Right = idxe(i,j,k) 
                Face(count)%qL => qBnd(:,4,:,:,Face(count)%Left ) 
                Face(count)%qR => qBnd(:,3,:,:,Face(count)%Right) 
                Face(count)%FL => FBnd(:,4,:,:,Face(count)%Left ) 
                Face(count)%FR => FBnd(:,3,:,:,Face(count)%Right) 
              ENDIF 
              Face(count)%nv = (/ 0., 1., 0. /) ! set face normal vector 
          ENDDO
         ENDDO
        ENDDO 
    ENDIF 
    ! z faces 
    IF(nDim>=3) THEN
        DO k = 1, KMAX+1
         DO j = 1, JMAX 
          DO i = 1, IMAX  
              count = count + 1 
              IF(k.EQ.1) THEN
                Face(count)%Left  = 0 
                Face(count)%Right = idxe(i,j,k)
                ALLOCATE( Face(count)%qL(nVar,nDOF(2),nDOF(3)) ) 
                ALLOCATE( Face(count)%FL(nVar,nDOF(2),nDOF(3)) ) 
                Face(count)%qR => qBnd(:,5,:,:,Face(count)%Right) 
                Face(count)%FR => FBnd(:,5,:,:,Face(count)%Right) 
              ELSEIF(k.EQ.KMAX+1) THEN 
                Face(count)%Left  = idxe(i,j,k-1) 
                Face(count)%Right = 0 
                Face(count)%qL => qBnd(:,6,:,:,Face(count)%Left ) 
                Face(count)%FL => FBnd(:,6,:,:,Face(count)%Left ) 
                ALLOCATE( Face(count)%qR(nVar,nDOF(2),nDOF(3)) ) 
                ALLOCATE( Face(count)%FR(nVar,nDOF(2),nDOF(3)) ) 
              ELSE              
                Face(count)%Left  = idxe(i,j,k-1) 
                Face(count)%Right = idxe(i,j,k) 
                Face(count)%qL => qBnd(:,6,:,:,Face(count)%Left ) 
                Face(count)%qR => qBnd(:,5,:,:,Face(count)%Right) 
                Face(count)%FL => FBnd(:,6,:,:,Face(count)%Left ) 
                Face(count)%FR => FBnd(:,5,:,:,Face(count)%Right) 
              ENDIF 
              Face(count)%nv = (/ 0., 0., 1. /) ! set face normal vector 
          ENDDO
         ENDDO
        ENDDO 
    ENDIF 
    !      
    DEALLOCATE( idxe ) 
    !
    ! We now need to define our basis functions. This is done by choosing a set of distinct 1D nodes in the interval [0,1] 
    ! The basis functions will be the Lagrange interpolation polynomials running through these nodes 
    CALL gauleg(0.,1.,xiGPN,wGPN,N+1)
    xin = xiGPN  ! WE make the following choice: the basis functions run through the Gauss-Legendre nodes (=> really orthogonal basis) 
    !
    ! Now, let us compute some of the important matrices in the ADER-DG method 
    ! 
    MM   = 0.   ! Element mass matrix 
    Kxi  = 0.   ! Element stiffness matrix 
    dudx = 0.   ! discrete derivative operator, which projects the derivatives onto the basis  
    DO iGP = 1, N+1 
        CALL BaseFunc1D(phi,phi_xi,xiGPN(iGP))
        DO k = 1, N+1
            DO l = 1, N+1
                ! i) Mass-matrix 
                MM(k,l) = MM(k,l) + wGPN(iGP)*phi(k)*phi(l) 
                ! ii) Stiffness-matrix 
                Kxi(k,l) = Kxi(k,l) + wGPN(iGP)*phi_xi(k)*phi(l)  
            ENDDO
        ENDDO        
    ENDDO    
    CALL MatrixInverse(N+1,MM,iMM) 
    dudx = MATMUL( iMM, TRANSPOSE(Kxi) ) 
    CALL BaseFunc1D(phi0,phi_xi,0.0) ! Compute the basis functions on the left 
    CALL BaseFunc1D(phi1,phi_xi,1.0) ! Compute the basis function on the right 
    ! The flux matrices are all possible combinations of left and right 
    DO k = 1, N+1
        DO l = 1, N+1
            FLm(k,l) = phi0(k)*phi1(l)   ! Left contribution to the left flux matrix    (m = left  of the interface)  
            FLp(k,l) = phi0(k)*phi0(l)   ! Right contribution to the left flux matrix   (p = right of the interface) 
            FRm(k,l) = phi1(k)*phi1(l)   ! Left contribution to the right flux matrix   (m = left  of the interface) 
            FRp(k,l) = phi1(k)*phi0(l)   ! Right contribution to the right flux matrix  (p = right of the interface) 
        ENDDO
    ENDDO        
    ! The time flux matrices for the ADER-DG predictor method are given by the principle of upwinding in time (causality principle) 
    F0 = phi0   ! upwinding in time = information comes from smaller times 
    F1 = FRm    ! upwinding in time = information comes from smaller times  
    K1 = F1 - Kxi 
    CALL MatrixInverse(N+1,K1,iK1)   
    FLcoeff = phi0  ! coefficients needed to extrapolate data onto the left  boundary 
    FRcoeff = phi1  ! coefficients needed to extrapolate data onto the right boundary 
    !
    ! For the fine output of each spatial degree of freedom onto a subgrid... 
    DO i = 1, N+1 
       subxi(i) = REAL(i-1)/REAL(N) 
    ENDDO
    cnt = 0 
    DO k = 1, N+1
     DO j = 1, N+1 
      DO i = 1, N+1  
         cnt = cnt + 1 
         CALL BaseFunc1D(phi_i,phi_xi,subxi(i))
         CALL BaseFunc1D(phi_j,phi_xi,subxi(j))
         CALL BaseFunc1D(phi_k,phi_xi,subxi(k))
         count = 0 
         DO kk = 1, nDOF(3) 
          DO jj = 1, nDOF(2)  
           DO ii = 1, nDOF(1) 
             count = count + 1 
             aux = (/ phi_i(ii), phi_j(jj), phi_k(kk) /) 
             SubOutputMatrix(count,cnt) = PRODUCT( aux(1:nDim) )                      
           ENDDO
          ENDDO
         ENDDO
         !
       ENDDO
      ENDDO
     ENDDO
     ! 
     ! subgrid triangulation for subgrid output 
     ! 
     ALLOCATE( idxn(N+1,N+1,N+1) )
     idxn = 0 
     c = 0 
     DO k = 1, N+1
        DO j = 1, N+1 
           DO i = 1, N+1 
              c = c + 1 
              allsubxi(:,c) = (/ REAL(i-1)/REAL(N), REAL(j-1)/REAL(N), REAL(k-1)/REAL(N) /)  
              idxn(i,j,k) = c
           ENDDO
        ENDDO
     ENDDO
     c = 0 
     DO k = 1, N 
        DO j = 1, N
           DO i = 1, N 
              c = c + 1 
              subtri(:,c) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i+1,j+1,k), idxn(i,j+1,k), idxn(i,j,k+1), idxn(i+1,j,k+1), idxn(i+1,j+1,k+1), idxn(i,j+1,k+1) /)         
           ENDDO
        ENDDO
     ENDDO
     DEALLOCATE( idxn )     
    !
    !
    ! Set the initial condition. Here, we assume a nodal basis. Otherwise, we would have to do formal L2 projection, 
    ! i.e. integration of the initial condition and multiplication with the inverse mass matrix. 
    !
    DO iElem = 1, nElem
        x0 = x(:,tri(1,iElem)) ! get the coordinate of the lower left node 
        DO k = 1, nDOF(3) 
         DO j = 1, nDOF(2)
          DO i = 1, nDOF(1) 
              xGP = x0 + (/ xiGPN(i), xiGPN(j), xiGPN(k) /)*dx(:) 
              CALL InitialField(u0,xGP) 
              uh(:,i,j,k,iElem) = u0 

          ENDDO
         ENDDO
        ENDDO
    ENDDO    
    !
    !CALL WriteDataTecplot
    CALL WriteDataGnuplot
    !
    CONTINUE 
    !
END SUBROUTINE ADERDGInit
    
    
SUBROUTINE InitialField(u0,xGP) 
    USE typesDef
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN ) :: xGP(d)             ! spatial position vector 
    REAL, INTENT(OUT) :: u0(nVar)           ! initial data vector in terms of conserved variables 
    ! Local variables 
    REAL :: VBase(nVar), ampl(nVar), sigma(d) 
    REAL :: V0(nVar) 
    REAL :: delta_rho, delta_T, delta_p, delta_vx, delta_vy, r
    REAL :: PI, epsilon
    ! 

    ! Gaussian perturbation 
!!$    sigma = (/ 0.05, 0.05, 0.05 /)       ! half-width
!!$    VBase(:) = (/ 1., 0., 0., 0., 1. /)  ! base-state 
!!$    ampl(:)  = 0.                        ! perturbation amplitude vector 
!!$    ampl(5)   = 1e-3                     ! 
!!$    V0(:) = VBase(:) + ampl(:)*EXP( -0.5*SUM(xGP(1:nDim)**2/sigma(1:nDim)**2) )    
    
!!$    PI = ACOS(-1.0)
!!$    epsilon = 5.0
!!$
!!$    r         = SQRT(( xGP(1) - 5.0 - time)**2 + ( xGP(2) - 5.0 - time)**2)
!!$    delta_T   = - epsilon**2*(EQN%gamma  -  1.0)/(8.0*EQN%gamma*PI**2)*EXP(1.0 - r**2)
!!$    delta_rho = (1.0 + delta_T)**(1.0/(EQN%gamma -  1.0)) - 1.0
!!$    delta_vx  = - ( xGP(2) - 5.0 - time)*epsilon/(2.0*PI)*EXP(0.5*(1.0 - r**2))
!!$    delta_vy  =   ( xGP(1) - 5.0 - time)*epsilon/(2.0*PI)*EXP(0.5*(1.0 - r**2))
!!$    delta_p   = (1.0 + delta_T)**(EQN%gamma/(EQN%gamma -  1.0)) - 1.0
!!$
!!$    V0(1) = 1.0 + delta_rho
!!$    V0(2) = 1.0 + delta_vx
!!$    V0(3) = 1.0 + delta_vy
!!$    V0(4) = 0.0
!!$    V0(5) = 1.0 + delta_p

    IF ( xGP(1).LT.0.5) then

    V0(1) = 1.0 
    V0(2) = 0.0 
    V0(3) = 0.0 
    V0(4) = 0.0
    V0(5) = 1.0 

    ELSE

    V0(1) = 0.1 
    V0(2) = 0.0 
    V0(3) = 0.0 
    V0(4) = 0.0
    V0(5) = 0.1 

    ENDIF

    ! A simple debug check for the computation of derivatives 
    !u0 = 0. 
    !u0(1) = 0.123 !*xGP(1) 
    CALL PDEPrim2Cons(u0,V0) 
    !
END SUBROUTINE InitialField
    
    
