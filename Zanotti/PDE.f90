SUBROUTINE PDEFlux(F,Q)
    USE typesDef, ONLY : nVar, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar)
    REAL, INTENT(OUT) :: F(nVar,d) 
    ! Local variables 
    REAL :: p, irho, lam, mu 
    !
    !
    ! 3D compressible Euler equations 
    !
    irho = 1.0/Q(1)
    p = (EQN%gamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)*irho )
    ! 
    F(1,1) = Q(2) 
    F(2,1) = irho*Q(2)*Q(2) + p 
    F(3,1) = irho*Q(2)*Q(3)
    F(4,1) = irho*Q(2)*Q(4)
    F(5,1) = irho*Q(2)*(Q(5)+p)  
    !
    F(1,2) = Q(3) 
    F(2,2) = irho*Q(3)*Q(2)  
    F(3,2) = irho*Q(3)*Q(3) + p 
    F(4,2) = irho*Q(3)*Q(4)
    F(5,2) = irho*Q(3)*(Q(5)+p)  
    ! 
    F(1,3) = Q(4) 
    F(2,3) = irho*Q(4)*Q(2)  
    F(3,3) = irho*Q(4)*Q(3)  
    F(4,3) = irho*Q(4)*Q(4) + p
    F(5,3) = irho*Q(4)*(Q(5)+p)  
    !
   
END SUBROUTINE PDEFlux
    
SUBROUTINE PDEEigenvalues(Lambda,Q,nv)
    USE typesDef, ONLY : nVar, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), nv(d) 
    REAL, INTENT(OUT) :: Lambda(nVar) 
    ! Local variables 
    REAL :: p, u, c, cp, cs, rho0, lam, mu  
    !
    u = ( Q(2)*nv(1) + Q(3)*nv(2) + Q(4)*nv(3) )/Q(1)       ! normal velocity 
    p = (EQN%gamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)/Q(1) )    ! fluid pressure 
    c = SQRT(EQN%gamma*p/Q(1))                              ! sound speed
    !
    Lambda = (/ u-c, u, u, u, u+c /)                        ! The eigenvalues of the Euler equations 
    !
    !
END SUBROUTINE PDEEigenvalues

SUBROUTINE PDECons2Prim(V,Q,iErr)
    USE typesDef, ONLY : nVar, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)     :: Q(nVar)     ! vector of conserved quantities 
    REAL, INTENT(OUT)    :: V(nVar)     ! primitive variables 
    INTEGER, INTENT(OUT) :: iErr        ! error flag 
    ! Local variables 
    REAL                 :: p 
    !    
    iErr = 0     
    !
    p = (EQN%gamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)/Q(1) )    ! fluid pressure 
    V(1) = Q(1)             ! fluid density 
    V(2:4) = Q(2:4)/Q(1)    ! fluid velocity 
    V(5)   = p              ! fluid pressure 
    ! 
    !
END SUBROUTINE PDECons2Prim 
        
SUBROUTINE PDEPrim2Cons(Q,V)
    USE typesDef, ONLY : nVar, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)      :: V(nVar)     ! primitive variables 
    REAL, INTENT(OUT)     :: Q(nVar)     ! vector of conserved quantities 
    ! Local variables 
    !   
    Q(1)   = V(1)           ! fluid density 
    Q(2:4) = V(1)*V(2:4)    ! momentum 
    Q(5)   = V(5)/(EQN%gamma-1) + 0.5*V(1)*SUM(V(2:4)**2)   ! total energy = internal energy + kinetic energy 
    !
END SUBROUTINE PDEPrim2Cons      
            
SUBROUTINE PDEVarName(Name) 
    USE typesDef, ONLY : nVar, d, EQN 
    IMPLICIT NONE       
    ! Argument list 
    CHARACTER(LEN=10), INTENT(OUT) :: Name(nVar)
    !
    Name(1) = 'rho' 
    Name(2) = 'u' 
    Name(3) = 'v'
    Name(4) = 'w' 
    Name(5) = 'p' 
    !
    !
END SUBROUTINE PDEVarName
    
    
    
    