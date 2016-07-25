!*****************************************************************************!
!* 
!*      $RCSfile: gauleg.f90,v $
!*
!*      $Revision: 4.2 $
!*
!*      $Date: 2006/04/07 09:24:50 $
!*
!*      $Author: iagmidu $
!*
!*****************************************************************************!

PURE SUBROUTINE gauleg(x1,x2,x,w,n)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  INTEGER     ::  n
  REAL        :: x1,x2,x(n),w(n)
  REAL        :: EPS
  INTEGER     :: i,j,m
  REAL        :: p1,p2,p3,pp,xl,xm,z,z1,Pi 
  !--------------------------------------------------------------------------
  INTENT(IN)  :: x1,x2,n
  INTENT(OUT) :: x,w
  !--------------------------------------------------------------------------
  PARAMETER (EPS=1E-15)
  !--------------------------------------------------------------------------

  m  = (n+1)/2
  xm = 0.5*(x2+x1)
  xl = 0.5*(x2-x1)
  Pi = ACOS(-1.0) 
  DO i=1,m
     z = COS(Pi*(i-.25)/(n+.5))
1    CONTINUE
     p1 = 1.
     p2 = 0.
     DO j = 1,n
        p3 = p2
        p2 = p1
        p1 = ((2.*j-1.)*z*p2-(j-1.)*p3)/j
     END DO
     pp = n*(z*p1-p2)/(z*z-1.)
     z1 = z
     z  = z1-p1/pp
     IF(ABS(z-z1).GT.EPS)GOTO 1
     x(i)    = xm-xl*z
     x(n+1-i)= xm+xl*z
     w(i)    = 2.*xl/((1.-z*z)*pp*pp)
     w(n+1-i)= w(i)
  END DO
  RETURN
END SUBROUTINE gauleg

!  (C) Copr. 1986-92 Numerical Recipes Software #

!
! Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde 
! matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
! integration and spectral methods. 
!
! Reference on LGL nodes and weights: 
!   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
!   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
!
! Written by Greg von Winckel - 04/17/2004
! Contact: gregvw@chtm.unm.edu
!
PURE SUBROUTINE gaulob(x1,x2,x,w,n1)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  ! Argument list  
  INTEGER     :: n1  
  REAL        :: x1,x2,x(n1),w(n1)
  ! Local variables 
  REAL        :: EPS
  INTEGER     :: i,k
  REAL        :: xold(n1), Pi  
  REAL        :: P(n1,n1) 
  INTEGER     :: n, iter, maxiter 
  !--------------------------------------------------------------------------
  INTENT(IN)  :: x1,x2,n1
  INTENT(OUT) :: x,w
  !--------------------------------------------------------------------------
  PARAMETER (EPS=1E-15)
  !--------------------------------------------------------------------------
  ! 
  n=N1-1;
  ! Use the Chebyshev-Gauss-Lobatto nodes as the first guess
  Pi = ACOS(-1.0) 
  DO i = 0, N 
   x(i+1) = cos(pi*REAL(i)/REAL(N)) 
  ENDDO 
  ! The Legendre Vandermonde Matrix
  P=0. 
  ! 
  ! Compute P_(N) using the recursion relation
  ! Compute its first and second derivatives and 
  ! update x using the Newton-Raphson method.
  ! 
  xold=2 
  !
  maxiter = 100000  
  DO iter = 1, maxiter 
    xold=x       
    P(:,1)=1
    P(:,2)=x   
    DO k = 2, N
        P(:,k+1)=( (2*k-1)*x*P(:,k)-(k-1)*P(:,k-1) )/REAL(k)
    ENDDO      
    x=xold-( x*P(:,N1)-P(:,N) )/( N1*P(:,N1) )             
    IF(MAXVAL(ABS(x-xold)).LE.eps) THEN
      EXIT
    ENDIF  
  ENDDO 
  w = 2.0/(N*N1*P(:,N1)**2)
  w = 0.5*w*(x2-x1) 
  x = x1 + 0.5*(x+1.)*(x2-x1)  
  ! 
END SUBROUTINE gaulob

