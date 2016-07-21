	MODULE M_PARS
	implicit none
	
!define parameters !!!!!!!!!!!!!!!!!!!!!!!
!m, q would be mass and q isn't used at the moment
!Rmin, Rout, left and right boundary location
!cfl, factor to specify ratio of dt/dx
!disip, value of the dissipation parameter
!tau, for a special way of handling the boundaries
!Nx, Nt: number of grid points and steps to be taken in time
!freq: how often to do output
!bc: choose between possible boundary condition

	real*8 :: m, q, Rmin, Rout, cfl, disip, tau
	integer :: Nx, Nt, freq, bc, derorder, disorder
	
	
	CONTAINS
	
	
	subroutine readpars
	implicit none
        namelist /pars_input/ Nx, Nt, freq, cfl, disip,derorder,disorder,&
    &                          m, q, Rmin, &
    &	                      Rout, bc, tau

!read params !!!!!!!!!!!!!!!!!!!!!!!


   open (unit = 10, file = "pars.in", status = "old" )
   read (unit = 10, nml = pars_input)
   close(unit = 10)
	
	end subroutine readpars
	
	
	end module M_pars	
