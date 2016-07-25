PROGRAM ADERDG3D 
  USE typesDef
  IMPLICIT NONE
  ! Local variables
  INTEGER :: i,j,k,iElem,iFace
  REAL, POINTER :: lQbndL(:,:,:),lFbndL(:,:,:),lQbndR(:,:,:),lFbndR(:,:,:)
  PRINT *, ' ----------------------------------------------- ' 
  PRINT *, '       A simple introduction to ADER-DG          ' 
  PRINT *, '          Written by and Copywright by           '
  PRINT *, '  Prof. Michael Dumbser and Dr. Olindo Zanotti   '
  PRINT *, '         University of Trento, Italy             ' 
  PRINT *, ' ----------------------------------------------- ' 

  ! We first need to compute the relevant matrices, set initial
  ! conditions and prepare some necessary stuff...  
  CALL ADERDGInit 
  CALL CPU_TIME(tCPU1) 
  ! Main loop in time 
  DO timestep = 1, NMAX
     IF(time >= tend) THEN
        EXIT
     ENDIF
     ! Compute the time step size according to the CFL condition 
     CALL CalcTimeStep  
     ! ADER predictor step 
     DO iElem  = 1, nElem
        CALL ADERSpaceTimePredictorNonlinear(qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem),qBnd(:,:,:,:,iElem),FBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem))  
     ENDDO
     ! Compute the element volume integral 
     DO iElem  = 1, nElem
        CALL ADERVolumeIntegral(duh(:,:,:,:,iElem),qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem))  
     ENDDO
     ! Set the boundary conditions (could depend on space and time)  
     CALL BoundaryConditions 
     ! Solve the Riemann problems and compute the surface integrals 
     DO iFace = 1, nFace
        CALL ADERRiemannSolver( Face(iFace)%qL,Face(iFace)%FL,Face(iFace)%qR,Face(iFace)%FR,Face(iFace)%nv )   
     ENDDO
     ! Compute the surface integrals of the test function multiplied with the numerical flux 
     DO iElem  = 1, nElem
        CALL ADERSurfaceIntegral(duh(:,:,:,:,iElem),FBnd(:,:,:,:,iElem))
     ENDDO
     ! Do the element update 
     DO iElem  = 1, nElem
        CALL ElementUpdate(uh(:,:,:,:,iElem),duh(:,:,:,:,iElem))
     ENDDO
     IF(MOD(timestep,10)==0) THEN
        PRINT *, ' n = ', timestep, ' t = ', time 
     ENDIF
     time = time + dt 

     CALL PlotNorm
     IF (  MOD ( timestep, 1) .EQ. 0) THEN  
        CALL WriteDataGnuplot
     ENDIF
  ENDDO
  CALL CPU_TIME(tCPU2)

  TEU = timestep*nElem 
  PRINT *, ' Total CPU time = ', tCPU2-tCPU1 
  PRINT *, ' Time / element update = ', (tCPU2-tCPU1)/TEU 
  PRINT *, ' Time / DOF update = ', (tCPU2-tCPU1)/TEU/PRODUCT(nDOF(1:nDim))  

  CALL WriteDataGnuplot
  CALL AnalyseError

  PRINT *, ' ----------------------------------------- ' 
  PRINT *, '  Program terminated!                      ' 
  PRINT *, ' ----------------------------------------- ' 


END PROGRAM ADERDG3D
    
