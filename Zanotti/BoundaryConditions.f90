SUBROUTINE BoundaryConditions 
  USE typesDef
  ! Local variables 
  REAL :: j,k,iFace
  REAL :: Qbc(nVar),Fbc(nVar,d),Vbc(nVar)
  !
  ! Fix boundary data  
  !Vbc = (/ 1., 0., 0., 0., 1. /)    ! primitive variables     
  !CALL PDEPrim2Cons(qBC,Vbc)        ! convert into conservative variables    
  !
  DO iFace = 1, nFace
     ! Here, we need to take care of the boundary conditions 
     ! For the moment, we use either simple extrapolation (copy from inside the domain) 
     ! or impose a constant value 
     IF(Face(iFace)%Left.EQ.0) THEN
        DO k = 1, nDOF(3)
           DO j = 1, nDOF(2) 
              Face(iFace)%qL(:,j,k) = Face(iFace)%qR(:,j,k)
              Face(iFace)%FL(:,j,k) = Face(iFace)%FR(:,j,k)
           ENDDO
        ENDDO
     ENDIF
     IF(Face(iFace)%Right.EQ.0) THEN 
        DO k = 1, nDOF(3)
           DO j = 1, nDOF(2) 
              Face(iFace)%qR(:,j,k) = Face(iFace)%qL(:,j,k)
              Face(iFace)%FR(:,j,k) = Face(iFace)%FL(:,j,k)
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  ! 
END SUBROUTINE BoundaryConditions
    
    
