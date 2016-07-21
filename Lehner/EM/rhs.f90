 	MODULE M_RHS
	implicit none
	
	contains
	
	subroutine rhs(u,dxu,dyU,dzU,source,dx,Nx,x)
	use m_pars
	implicit none
	real*8, dimension(:,:,:,:), intent(in) :: u, dxu,dyu, dzU
	real*8, dimension(:,:,:,:), intent(inout) :: source
	real*8, dimension(:,:,:), intent(in) :: x
	real*8 :: dx, v
	integer :: i, Nx
	logical :: ltrace


!for debug
	ltrace=.true.

	
	source = 0.0
		
!electric field eqns + damp
	source(iEx,:,:,:) =  dyU(iBz,:,:,:) - dzU(iBy,:,:,:) - dxU(iff,:,:,:)*k_r 

	source(iEy,:,:,:) = -dxU(iBz,:,:,:) + dzU(iBx,:,:,:) - dyU(iff,:,:,:)*k_r 
	
	source(iEz,:,:,:) =  dxU(iBy,:,:,:) - dyU(iBx,:,:,:) - dzU(iff,:,:,:)*k_r

	source(iff,:,:,:) = -dxU(iEx,:,:,:)-dyU(iEy,:,:,:)-dzU(iEz,:,:,:) - k_r*u(iff,:,:,:)


!magnetic field eqns + damp
	source(iBx,:,:,:) = -(dyU(iEz,:,:,:) - dzU(iEy,:,:,:)) - dxU(igg,:,:,:) *k_r

	source(iBy,:,:,:) = -(-dxU(iEz,:,:,:) + dzU(iEx,:,:,:)) - dyU(igg,:,:,:) *k_r
	
	source(iBz,:,:,:) = -(dxU(iEy,:,:,:) - dyU(iEx,:,:,:)) - dzU(igg,:,:,:)*k_r

	source(igg,:,:,:) = -dxU(iBx,:,:,:)-dyU(iBy,:,:,:)-dzU(iBz,:,:,:) - k_r*u(igg,:,:,:)

			
	end subroutine rhs
	
	end module M_RHS
