	MODULE M_PARS
	implicit none
	
!define parameters !!!!!!!!!!!!!!!!!!!!!!!
	real*8 ::  cfl, disip, sL, sR, c,s,a,k_r,k_e, gamma, cl,valuefloor
	integer :: N, Nt, freq, bc, derorder,dissorder,STMETH

	integer :: iEx,iEy,iEz,iff, iBx,iBy,iBz,igg	
	
	CONTAINS
	
	
	subroutine readpars
	implicit none
        namelist /pars_input/ N, Nt, freq, cfl, disip, &
	&		sL, sR, derorder,dissorder, c,s,a, STMETH, &
	&		k_r, k_e, gamma, cl,valuefloor

!read params !!!!!!!!!!!!!!!!!!!!!!!


   open (unit = 10, file = "pars.in", status = "old" )
   read (unit = 10, nml = pars_input)
   close(unit = 10)


!electric field and damping field ff
	iEx = 1
	iEy = 2
	iEz = 3
	iff = 4


!magnetic field and damping field gg
	iBx = 5
	iBy = 6
	iBz = 7
	igg = 8


	
	end subroutine readpars
	
	
	end module M_pars	
