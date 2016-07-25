!
! Small and simple Gauss package
! written by Michael Dumbser
! 27.05.2011 
!
SUBROUTINE LinSolve(N,A,b,x)
  IMPLICIT NONE
  INTEGER       :: N
  REAL          :: A(N,N), b(N), x(N)  
  !
  INTEGER       :: i,j,flag,ml(1) 
  REAL          :: temp(N+1),piv 
  REAL          :: C(N+1,N)  
  !
  C(1:N,:)     = TRANSPOSE(A)
  C(N+1,:)     = b 
  !    
  ! Forward elimination with column pivoting 
  ! 
  DO i = 1, N
     ! If pivot element is zero, then swap rows 
     ml = MAXLOC(ABS(C(i,i:N))) 
     j = i - 1 + ml(1) 
     temp   = C(:,j) 
     C(:,j) = C(:,i)
     C(:,i) = temp      
     IF(C(i,i).EQ.0.) THEN
        PRINT *, 'ERROR. Matrix is singular!'
        DO j = 1, N
           PRINT *, A(j,:) 
        ENDDO
        STOP
     ENDIF
     piv    = 1./C(i,i)
     C(:,i) = C(:,i)*piv 
     DO j = i+1, N 
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  ! Back substitution
  !
  DO i = N,1,-1   
     DO j = i-1,1,-1
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  x = C(N+1,:)  
  !
END SUBROUTINE LinSolve

SUBROUTINE LinSolve2(N,M,A,b,x,iErr)
   IMPLICIT NONE
   INTEGER       :: N,M
   INTEGER       :: iErr 
   REAL          :: A(N,N), b(N,M), x(N,M)
   !
   INTEGER       :: i,j,ml(1)
   REAL          :: temp(N+M),piv
   REAL          :: C(N+M,N)
   !
   iErr = 0 
   ! 
   C(1:N,:)      = TRANSPOSE(A)
   C(N+1:N+M,:)  = TRANSPOSE(b)
   !
   ! Forward elimination with column pivoting
   !
   DO i = 1, N
     ! If pivot element is zero, then swap rows
     ml = MAXLOC(ABS(C(i,i:N)))
     j = i - 1 + ml(1)
     temp   = C(:,j)
     C(:,j) = C(:,i)
     C(:,i) = temp
     IF(C(i,i).EQ.0.) THEN
        PRINT *, 'ERROR. Matrix is singular!'
        iErr = 1 
        RETURN 
     ENDIF
     piv    = 1./C(i,i)
     C(:,i) = C(:,i)*piv
     DO j = i+1, N
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
   ENDDO
   !
   ! Back substitution
   !
   DO i = N,1,-1
    DO j = i-1,1,-1
      C(:,j) = C(:,j) - C(i,j)*C(:,i)
    ENDDO
   ENDDO
   !
   x = TRANSPOSE(C(N+1:N+M,:))
   !
END SUBROUTINE LinSolve2    
    
SUBROUTINE MatrixInverse(N,A,iA)
  IMPLICIT NONE
  INTEGER       :: N
  REAL          :: A(N,N), iA(N,N)
  !
  INTEGER       :: i,j,flag,ml(1) 
  REAL          :: piv
  REAL          :: temp(2*N)
  REAL          :: C(2*N,N)  
  !
  C(1:N,:)     = TRANSPOSE(A)
  C(N+1:2*N,:) = 0. 
  DO i = 1, N
     C(N+i,i) = 1.
  ENDDO
  !    
  ! Forward elimination and row swapping (if necessary)
  ! 
  DO i = 1, N
     ! If pivot element is zero, then swap rows 
     ml = MAXLOC(ABS(C(i,i:N))) 
     j = i - 1 + ml(1) 
     temp   = C(:,j) 
     C(:,j) = C(:,i)
     C(:,i) = temp      
     IF(C(i,i).EQ.0.) THEN
        PRINT *, 'ERROR. Matrix is singular!'
        DO j = 1, N
           PRINT *, A(j,:) 
        ENDDO
        STOP
     ENDIF
     piv    = 1./C(i,i)
     C(:,i) = C(:,i)*piv 
     DO j = i+1, N 
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  ! Back substitution
  !
  DO i = N,1,-1   
     DO j = i-1,1,-1
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  iA = TRANSPOSE( C(N+1:2*N,:) ) 
  !
END SUBROUTINE MatrixInverse
