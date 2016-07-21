	MODULE M_BOUNDARY
	use m_pars
	
	implicit none
	
	contains
	
	
	subroutine boundary(u,du,uboun,dx)
	real*8, dimension(:,:), intent(in) :: u, du
	real*8, dimension(:,:), intent(inout) :: uboun
	real*8 dx
	
!Olsson
	if(bc.eq.0) then
	uboun(2,nx) = 0.0
	uboun(3,nx) = 0.0

!carpenter
	else if(bc.eq.1) then
	uboun(:,nx) = uboun(:,nx) - tau * u(:,nx)/dx
	
!periodic			
	else if(bc.eq.-1) then
	uboun(:,nx) = (uboun(:,1)+uboun(:,nx))*0.5
	uboun(:,1) = uboun(:,nx) 


!fully reflecting
	else if(bc.eq.2) then
	uboun(2,1) = 0.0
	uboun(2,nx) = 0.0
	uboun(3,1) = du(2,1)+du(3,1)
 	uboun(3,nx) = du(2,nx)-du(3,nx)
	end if

	
	end subroutine boundary
	
	
	end module M_BOUNDARY
	
	
	
