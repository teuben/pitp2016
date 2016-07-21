	MODULE M_RHS
	implicit none
	
	contains
	
	subroutine rhs(u,du,rk,ru,conduc,dx,Nx, eb,x)
	use m_pars, only : cfl, c
	
	implicit none
	real*8, dimension(:,:), intent(in) :: u, du
	real*8, dimension(:,:), intent(inout) :: rk,ru
	real*8, dimension(:) :: x, conduc
	real*8 :: dx,xx,vv, bb,pen
	integer :: i, Nx, eb, TT, T1, T2
	
	rk = 0.0
	ru = 0.0

        ! u(1,i) = V    ,   u(2, i) = W

	do i=1, Nx	
   	  rk(1,i) =  du(2,i)
   	  rk(2,i) =  du(1,i) 
	  ru(1,i) = -conduc(i)*(u(1,i) - c * u(2,i))
	  ru(2,i) = 0.0d0
	end do


	end subroutine rhs
	
	end module M_RHS
