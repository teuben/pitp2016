	module m_initial
	implicit none
	
	contains
	
	subroutine initial(uold,x, Nx)
	implicit none
	real*8, dimension(:,:) :: uold
	real*8, dimension(:) :: x
	integer Nx, i
	
	do i=1, Nx
	uold(:,i) = 0.0
	end do
	
!simple initial data take out
	do i=1,NX
	  if(x(i).ge.7.and.x(i).le.9.) then
	    uold(1,i) = (x(i)-7.)**4 * (x(i)-9.)**4 
	    uold(2,i) = 4.*((x(i)-7.)**3*(x(i)-9.)**4 &
     &                     +(x(i)-7.)**4 * (x(i)-9.)**3)
	    uold(3,i) = uold(2,i)
	  end if
	end do
	
	end subroutine initial
	
	end module m_initial
