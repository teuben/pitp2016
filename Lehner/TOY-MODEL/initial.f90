	module m_initial
	implicit none
	
	contains
	
	subroutine initial(uold,x, Nx)
	use m_pars, only : c
	implicit none
	real*8, dimension(:,:) :: uold
	real*8, dimension(:) :: x
	integer Nx, i
	real*8 bb, norma

	bb = 1.5
	norma = 2.0*18.966528
	
	do i=1, Nx
	  uold(:,i) = 0.0
	end do
	
!stupid initial data take out
	do i=1,NX
         if(x(i).ge.0.5.and.x(i).le.0.6) then
          uold(2,i) = 4.*(x(i)-0.5d0)**3*(x(i)-0.6d0)**3/0.05**8 &
                       *((x(i)-0.5d0)+(x(i)-0.6d0)) 
	  uold(2,i) = uold(2,i)/(norma*c)	       
          uold(1,i) = c * uold(2,i)
         end if
	end do
	
	end subroutine initial
	
	end module m_initial
