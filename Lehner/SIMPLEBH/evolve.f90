	module M_EVOLVE
	implicit none
	
	contains
	
	
	subroutine evolve(unew,uold,dx,dt,x)
	use m_derivs
	use m_boundary
	use m_rhs
	use m_pars, only: disip
	implicit none
	
	real*8, dimension(:,:), intent(inout):: unew
	real*8, dimension(:,:), intent(in):: uold
	real*8, dimension(:) :: x
	real*8 dx, dt
	
!keep in some variables the rhs
	real*8, dimension(:,:), allocatable :: du,urk1,urk2,urk3,diss
	
	integer i
	
!allocate memmory
	allocate(du(9,Nx),urk1(9,Nx),urk2(9,Nx),urk3(9,Nx),diss(9,nx))
	
! We'll do the rk 3rd stuff. The procedure is basically
! the same 3 times, (i) evaluate derivs; (ii) evaluate rhs;
! (iii) get the new value
	
!!!!!!!!!!!!!!!!!!! FIRST STEP	

	
!calculate derivs	
	do i = 1, 3
	  call derivs(uold(i,:),du(i,:),dx,Nx)
	  call dissip(uold(i,:),diss(i,:),dx,Nx)
	end do

	
	call rhs(uold,du,urk1,dx,Nx,x)

	call boundary(uold,du,urk1,dx)
	
	unew = uold + dt * 0.5 * (urk1+disip*diss)
	
!!!!!!!!!!!!!!!!!!! SECOND STEP	
	
!calculate derivs	
	do i = 1, 3
	  call derivs(unew(i,:),du(i,:),dx,Nx)
	  call dissip(unew(i,:),diss(i,:),dx,Nx)	  
	end do
	
	call rhs(unew,du,urk2,dx,Nx,x)
	
	call boundary(unew,du,urk2,dx)
	
	unew = uold + dt * 0.75 * (urk2+disip*diss)

!!!!!!!!!!!!!!!!!!! THIRD STEP	
	
!calculate derivs	
	do i = 1, 3
	  call derivs(unew(i,:),du(i,:),dx,Nx)
	  call dissip(unew(i,:),diss(i,:),dx,Nx)	  	  
	end do
	  
	call rhs(unew,du,urk3,dx,Nx,x)
	urk3 = urk3 + disip * diss
	
	call boundary(unew,du,urk3,dx)
	
	unew = uold + dt/9. * (2.*urk1+3.*urk2+4.*urk3)		
	
!deallocate memmory

	deallocate(du,urk3,urk2,urk1)

	
	end subroutine evolve
	
	
	end module M_EVOLVE
