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

! 1 = phi
! 2 = g==phi,t
! 3 = f== phi,x

	do i=1, Nx

	source(:,i) = 0.0	
	
	end do	
!exmaple here take out later!
	do i =1, Nx

	source(1,i) = u(2,i)
	source(2,i) = du(3,i)
	source(3,i) = du(2,i)

	end do

	end subroutine rhs
	
	end module M_RHS
