	MODULE M_INITIAL
	implicit none
	
	CONTAINS
	
	subroutine initial(Uold,dx,dy,dz,x,y,z,Nv,N)
	use m_pars, only :  iEx,iEy,iEz,iff,iBx,iBy,iBz,igg
	real*8, dimension(:,:,:,:) :: Uold
	real*8, dimension(:,:,:) :: x,y,z
	real*8 :: dx, dy, dz, r2
	integer :: Nv, N, i,j,k
	
	real*8 :: ppi

	ppi = 3.141592653589793

	uold = 0.0d0

	do i=1,N
	  do j=1,N
	   do k=1,N
	     r2 = sqrt(x(i,j,k)**2+y(i,j,k)**2)
	    if(r2.le.5.) then
	      uold(iBz,i,j,k) = (r2-5.0d0)**4 + 1e-5*sin(z(i,j,k))
	    end if
	  end do
	 end do
 	end do
	
	end subroutine initial
	
	
	END MODULE M_INITIAL
