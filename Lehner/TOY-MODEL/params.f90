	MODULE M_PARS
	implicit none
	
!define parameters !!!!!!!!!!!!!!!!!!!!!!!
	real*8 :: m, q, Rout, cfl, disip, cond, tau, time,loc, c
	integer :: dim, Nx, Nt, freq, bc, pos, derorder, runge_kutta
	
	
	CONTAINS
	
	
	subroutine readpars
	implicit none
        namelist /pars_input/ dim, runge_kutta, Nx, Nt, freq, cfl, disip, cond, c,  m, q, Rout, bc, tau, derorder

!read params !!!!!!!!!!!!!!!!!!!!!!!


   open (unit = 10, file = "pars.in", status = "old" )
   read (unit = 10, nml = pars_input)
   close(unit = 10)
	
	end subroutine readpars
	
	
	end module M_pars	
