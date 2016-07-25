SUBROUTINE ADERSurfaceIntegral(lduh,lFbnd)
    USE typesDef
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)    :: lFbnd(nVar,6,nDOF(2),nDOF(3))            ! nonlinear flux tensor in each space-time DOF 
    REAL, INTENT(INOUT) :: lduh(nVar,nDOF(1),nDOF(2),nDOF(3))       ! spatial degrees of freedom 
    ! Local variables 
    INTEGER           :: i,j,k,l,iVar 
    REAL              :: aux(d) 
    !
    ! Now multiply the numerical fluxes on the surfaces with the test functions and compute the surface integrals 
    ! 
    ! x faces
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2) 
            aux = (/ 1., wGPN(j), wGPN(k) /)
            DO iVar = 1, nVar 
                lduh(iVar,:,j,k) = lduh(iVar,:,j,k) - PRODUCT(aux(1:nDim))/dx(1)*( lFbnd(iVar,2,j,k)*FRCoeff - lFbnd(iVar,1,j,k)*FLCoeff )      ! left flux minus right flux 
            ENDDO                                                             
        ENDDO
    ENDDO 
    IF(nDim>=2) THEN
        ! y faces
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
                aux = (/ 1., wGPN(i), wGPN(k) /) 
                DO iVar = 1, nVar 
                    lduh(iVar,i,:,k) = lduh(iVar,i,:,k) - PRODUCT(aux(1:nDim))/dx(2)*( lFbnd(iVar,4,i,k)*FRCoeff - lFbnd(iVar,3,i,k)*FLCoeff )  ! left flux minus right flux  
                ENDDO                                                             
            ENDDO
        ENDDO 
    ENDIF
    IF(nDim>=3) THEN
        ! z faces
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1) 
                aux = (/ 1., wGPN(i), wGPN(j) /) 
                DO iVar = 1, nVar 
                    lduh(iVar,i,j,:) = lduh(iVar,i,j,:) - PRODUCT(aux(1:nDim))/dx(3)*( lFbnd(iVar,6,i,j)*FRCoeff - lFbnd(iVar,5,i,j)*FLCoeff )  ! left flux minus right flux  
                ENDDO                                                             
            ENDDO
        ENDDO 
    ENDIF
    !
END SUBROUTINE ADERSurfaceIntegral 
    
    