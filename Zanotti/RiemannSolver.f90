SUBROUTINE ADERRiemannSolver(lQbndL,lFbndL,lQbndR,lFbndR,nv)
  USE typesDef
  IMPLICIT NONE 
  ! Argument list 
  REAL, INTENT(IN)     :: lQbndL(nVar,nDOF(2),nDOF(3))              ! left space-time degrees of freedom 
  REAL, INTENT(IN)     :: lQbndR(nVar,nDOF(2),nDOF(3))              ! right space-time degrees of freedom 
  REAL, INTENT(IN)     :: nv(d)                                     ! normal vector 
  REAL, INTENT(INOUT)  :: lFbndL(nVar,nDOF(2),nDOF(3))              ! left flux 
  REAL, INTENT(INOUT)  :: lFbndR(nVar,nDOF(2),nDOF(3))              ! right flux 
  ! Local variables 
  INTEGER           :: i,j,k,l
  REAL              :: aux(d), QavL(nVar), QavR(nVar), smax
  REAL              :: LL(nVar), LR(nVar) 
  !
  ! Compute the average states from the left and the right, which we need to compute the numerical dissipation 
  QavL = 0. 
  QavR = 0. 
  DO k = 1, nDOF(3)
     DO j = 1, nDOF(2)
        aux = (/ 1., wGPN(j), wGPN(k) /) 
        QavL = QavL + PRODUCT(aux(1:nDim))*lQbndL(:,j,k) 
        QavR = QavR + PRODUCT(aux(1:nDim))*lQbndR(:,j,k) 
     ENDDO
  ENDDO
  !
  ! Here, we implement a very simple Rusanov scheme with scalar dissipation (smax*Id). 
  !
  CALL PDEEigenvalues(LL,QavL,nv) 
  CALL PDEEigenvalues(LR,QavR,nv) 
  smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) ) 
  !
  ! We now compute the numerical flux. Note that the scheme is at the moment written in 
  ! CONSERVATION FORM => no fluctuations, but real fluxes. 
  ! Later, this will be converted into the left and right fluctuations. 
  !
  DO k = 1, nDOF(3)
     DO j = 1, nDOF(2)
        lFbndL(:,j,k) = 0.5*( lFbndR(:,j,k) + lFbndL(:,j,k) ) - 0.5*smax*( lQbndR(:,j,k) - lQbndL(:,j,k) ) 
        lFbndR(:,j,k) = lFbndL(:,j,k) 
     ENDDO
  ENDDO
  !
END SUBROUTINE ADERRiemannSolver
    
    
