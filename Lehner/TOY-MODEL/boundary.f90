	MODULE M_BOUNDARY
	use m_pars
	
	implicit none
	
	contains
	
	
	subroutine boundary(u,du,uboun,dx)
	real*8, dimension(:,:), intent(in) :: u, du
	real*8, dimension(:,:), intent(inout) :: uboun
	real*8 dx, gs, factor
	
!Olsson
	if(bc.eq.0) then
!	uboun(1,nx) = 0./(1.+rout) 
	uboun(2,nx) = 0/(1.+rout)**2
	uboun(3,nx) = 0./(1.+rout)**2

!carpenter
	else if(bc.eq.1) then
        if(derorder.eq.1) then
        factor = 2.
        else if(derorder.eq.2) then
        factor = 48./17.
        else if (derorder.eq.3) then
        factor = 43200./13649.
        else if (derorder.ge.4) then
        factor =  5080320./1498139.
        end if


	gs = 0.1 * sin(time)
!	uboun(1,nx) = uboun(1,nx) - factor * (u(1,nx)-1./(1.+rout))/dx
!        uboun(2,nx) = uboun(2,nx) - factor * (u(2,nx)+1.5/(1.+rout)**2)/dx
!        uboun(3,nx) = uboun(3,nx) - factor * (u(3,nx)+1./(1.+rout)**2)/dx
	
!periodic			
	else if(bc.eq.-1) then
   	 uboun(:,nx) = (uboun(:,1)+uboun(:,nx))*0.5
	 uboun(:,1) = uboun(:,nx) 

	else if(bc.eq.-2) then
   	 uboun(:,nx) = uboun(:,2)
	 uboun(:,1)  = uboun(:,nx-1) 
	end if
	
	end subroutine boundary
	
	
	end module M_BOUNDARY
	
	
	
