SUBROUTINE ADERVolumeIntegral(lduh,lqhi,lFhi)
    USE typesDef
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)  :: lqhi(nVar,nDOF(1),nDOF(2),nDOF(3))      ! space-time degrees of freedom 
    REAL, INTENT(IN)  :: lFhi(nVar,d,nDOF(1),nDOF(2),nDOF(3))    ! nonlinear flux tensor in each space-time DOF 
    REAL, INTENT(OUT) :: lduh(nVar,nDOF(1),nDOF(2),nDOF(3))      ! spatial degrees of freedom 
    ! Local variables 
    INTEGER           :: i,j,k,l
    REAL              :: aux(d) 
    ! 
    ! Initialize the update DOF 
    lduh = 0. 
    ! x - direction 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2) 
            aux = (/ 1., wGPN(j), wGPN(k) /) 
            lduh(:,:,j,k) = lduh(:,:,j,k) + MATMUL( lFhi(:,1,:,j,k), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(1) 
        ENDDO
    ENDDO
    ! y - direction (not needed in 1D) 
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
                aux = (/ 1., wGPN(i), wGPN(k) /) 
                lduh(:,i,:,k) = lduh(:,i,:,k) + MATMUL( lFhi(:,2,i,:,k), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(2)
            ENDDO
        ENDDO
    ENDIF 
    ! z - direction (not needed in 1D and 2D) 
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                aux = (/ 1., wGPN(i), wGPN(j) /)  
                lduh(:,i,j,:) = lduh(:,i,j,:) + MATMUL( lFhi(:,3,i,j,:), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(3)  
            ENDDO
        ENDDO
    ENDIF 
    !
    CONTINUE
    !
END SUBROUTINE ADERVolumeIntegral 
    
    
