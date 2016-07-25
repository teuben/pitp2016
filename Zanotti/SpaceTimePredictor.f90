SUBROUTINE ADERSpaceTimePredictorNonlinear(lqhi,lFhi,lQbnd,lFbnd,luh)
  USE typesDef
  IMPLICIT NONE
  ! Argument list
  REAL, INTENT(IN)  :: luh(nVar,nDOF(1),nDOF(2),nDOF(3))              ! spatial degrees of freedom
  REAL, INTENT(OUT) :: lqhi(nVar,nDOF(1),nDOF(2),nDOF(3))             ! time-averaged space-time degrees of freedom
  REAL, INTENT(OUT) :: lFhi(nVar,d,nDOF(1),nDOF(2),nDOF(3))           ! time-averaged nonlinear flux tensor in each space-time DOF
  REAL, INTENT(OUT) :: lqbnd(nVar,6,nDOF(2),nDOF(3))                  ! time-averaged space-time degrees of freedom
  REAL, INTENT(OUT) :: lFbnd(nVar,6,nDOF(2),nDOF(3))                  ! time-averaged nonlinear flux tensor in each space-time DOF
  ! Local variables
  INTEGER :: i,j,k,l,iVar,iDim, iter
  REAL    :: rhs0(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))               ! contribution of the initial condition to the known right hand side
  REAL    :: rhs(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! known right hand side
  REAL    :: lqh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! space-time degrees of freedom
  REAL    :: lFh(nVar,d,nDOF(1),nDOF(2),nDOF(3),nDOF(0))              ! nonlinear flux tensor in each space-time DOF
  REAL    :: aux(d), w                                                ! auxiliary variables
  REAL    :: lqhold(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))             ! old space-time degrees of freedom
  REAL    :: lqt(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! time derivative qt of q
  REAL    :: res                                                      ! residual
  REAL, PARAMETER :: tol = 1e-7                                      ! tolerance
  INTEGER, PARAMETER :: MAXNEWTON = 50
  !
  DO k = 1, nDOF(3)
     DO j = 1, nDOF(2)
        DO i = 1, nDOF(1)
           ! Trivial initial guess (can be significantly improved)
           DO iVar = 1, nVar
              lqh(iVar,i,j,k,:) = luh(iVar,i,j,k)
           ENDDO
           ! Compute the contribution of the initial condition uh to the time update. I prefer to compute it once
           ! and store it in rhs0, but if you think it is faster, you can also recompute this contribution
           ! inside the Picard loop (DO iter = 1, N+1)
           aux = (/ wGPN(i), wGPN(j), wGPN(k) /)
           DO iVar = 1, nVar
              rhs0(iVar,i,j,k,:) = PRODUCT(aux(1:nDim))*F0(:)*luh(iVar,i,j,k)
           ENDDO
           !
        ENDDO
     ENDDO
  ENDDO
  !
  ! Discrete Picard iterations. 
  DO iter = 1, MAXNEWTON
     ! save old space-time DOF
     lqhold = lqh
     DO l = 1, nDOF(0) ! loop over DOF in time
        ! Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
        DO k = 1, nDOF(3)
           DO j = 1, nDOF(2)
              DO i = 1, nDOF(1)
                 CALL PDEFlux(lFh(:,:,i,j,k,l),lqh(:,i,j,k,l))
              ENDDO
           ENDDO
        ENDDO
        ! Compute the "derivatives" (contributions of the stiffness matrix)
        ! x direction (independent from the y and z derivatives)
        DO k = 1, nDOF(3)
           DO j = 1, nDOF(2)
              aux = (/ wGPN(l), wGPN(j), wGPN(k) /)
              rhs(:,:,j,k,l) = rhs0(:,:,j,k,l) - PRODUCT(aux(1:nDim))*dt/dx(1)*MATMUL( lFh(:,1,:,j,k,l), Kxi )
           ENDDO
        ENDDO
        ! y direction (independent from the x and z derivatives) - should not be used for 1D
        IF(nDim>=2) THEN
           DO k = 1, nDOF(3)
              DO i = 1, nDOF(1)
                 aux = (/ wGPN(l), wGPN(i), wGPN(k) /)
                 rhs(:,i,:,k,l) = rhs(:,i,:,k,l) - PRODUCT(aux(1:nDim))*dt/dx(2)*MATMUL( lFh(:,2,i,:,k,l), Kxi )
              ENDDO
           ENDDO
        ENDIF
        ! z direction (independent from the x and y derivatives) - should not be used for 1D and 2D
        IF(nDim>=3) THEN
           DO j = 1, nDOF(2)
              DO i = 1, nDOF(1)
                 aux = (/ wGPN(l), wGPN(i), wGPN(j) /)
                 rhs(:,i,j,:,l) = rhs(:,i,j,:,l) - PRODUCT(aux(1:nDim))*dt/dx(3)*MATMUL( lFh(:,3,i,j,:,l), Kxi )
              ENDDO
           ENDDO
        ENDIF
        !
     ENDDO ! end loop over time DOF
     !
     ! Multiply with (K1)^(-1) to get the discrete time integral of the discrete Picard iteration
     !
     DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
           DO i = 1, nDOF(1)
              aux = (/ wGPN(i), wGPN(j), wGPN(k) /)
              lqh(:,i,j,k,:) = 1./(PRODUCT(aux(1:nDim)))*MATMUL( rhs(:,i,j,k,:), TRANSPOSE(iK1) )
              !lqt(:,i,j,k,:) = 1.0/dt*MATMUL( lqh(:,i,j,k,:), TRANSPOSE(dudx) )         ! currently used only for debugging purposes, to check if derivatives are correctly computed
           ENDDO
        ENDDO
     ENDDO
     !
     ! We can stop the iterations if a certain tolerance has been reached. 
     !
     res = SQRT(SUM((lqh-lqhold)**2))
     IF(res.LT.tol) THEN
        EXIT
     ENDIF
     !
  ENDDO
  !
  ! Immediately compute the time-averaged space-time polynomials
  !
  DO k = 1, nDOF(3)
     DO j = 1, nDOF(2)
        DO i = 1, nDOF(1)
           lqhi(:,i,j,k) = MATMUL( lqh(:,i,j,k,:), wGPN )
           DO iDim = 1, nDim
              lFhi(:,iDim,i,j,k) = MATMUL( lFh(:,iDim,i,j,k,:), wGPN )
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  ! Compute the bounday-extrapolated values for Q and F*n
  !
  lQbnd = 0.
  lFbnd = 0.
  ! x-direction: face 1 (left) and face 2 (right)
  DO k = 1, nDOF(3)
     DO j = 1, nDOF(2)
        lQbnd(:,1,j,k) = MATMUL( lqhi(:,:,j,k),   FLCoeff )   ! left
        lQbnd(:,2,j,k) = MATMUL( lqhi(:,:,j,k),   FRCoeff )   ! right
        lFbnd(:,1,j,k) = MATMUL( lFhi(:,1,:,j,k), FLCoeff )   ! left
        lFbnd(:,2,j,k) = MATMUL( lFhi(:,1,:,j,k), FRCoeff )   ! right
     ENDDO
  ENDDO
  ! y-direction: face 3 (left) and face 4 (right)
  IF(nDim>=2) THEN
     DO k = 1, nDOF(3)
        DO i = 1, nDOF(1)
           lQbnd(:,3,i,k) = MATMUL( lqhi(:,i,:,k),   FLCoeff )   ! left
           lQbnd(:,4,i,k) = MATMUL( lqhi(:,i,:,k),   FRCoeff )   ! right
           lFbnd(:,3,i,k) = MATMUL( lFhi(:,2,i,:,k), FLCoeff )   ! left
           lFbnd(:,4,i,k) = MATMUL( lFhi(:,2,i,:,k), FRCoeff )   ! right
        ENDDO
     ENDDO
  ENDIF
  ! z-direction: face 5 (left) and face 6 (right)
  IF(nDim>=3) THEN
     DO j = 1, nDOF(2)
        DO i = 1, nDOF(1)
           lQbnd(:,5,i,j) = MATMUL( lqhi(:,i,j,:),   FLCoeff )   ! left
           lQbnd(:,6,i,j) = MATMUL( lqhi(:,i,j,:),   FRCoeff )   ! right
           lFbnd(:,5,i,j) = MATMUL( lFhi(:,3,i,j,:), FLCoeff )   ! left
           lFbnd(:,6,i,j) = MATMUL( lFhi(:,3,i,j,:), FRCoeff )   ! right
        ENDDO
     ENDDO
  ENDIF
  !
  CONTINUE
  !
END SUBROUTINE ADERSpaceTimePredictorNonlinear








