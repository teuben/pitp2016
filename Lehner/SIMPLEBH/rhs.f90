	MODULE M_RHS
	implicit none
	
	contains
	
	subroutine rhs(u,du,source,dx,Nx,x)
	implicit none
	real*8, dimension(:,:), intent(in) :: u, du
	real*8, dimension(:,:), intent(inout) :: source
	real*8, dimension(:) :: x
	real*8 :: dx
	integer :: i, Nx

	real*8 :: M, l, r, ppi, phi, drppi, drphi, g, drg

	M=1.
	l=2.

	do i=1, Nx

	source(:,i) = 0.0	
	
	end do	
!exmaple here take out later!
	do i =1, Nx
	r = x(i)

	ppi = u(3,i)
	phi = u(1,i)
	g = u(2,i)
        drppi = du(3,i)
        drphi = du(1,i)
        drg = du(2,i)

	source(3,i) = 2.*(r-M)*g &
      &               + r*(r-2.*M)*drg &
      &               +2.*M*ppi &
      &               +4.*r*M*drppi &
      &               -l*(l+1.)*phi
	source(3,i) = source(3,i)/(r**2+2.*M*r)
	source(1,i) = ppi
	source(2,i) = drppi


 

	end do

	end subroutine rhs
	
	end module M_RHS
