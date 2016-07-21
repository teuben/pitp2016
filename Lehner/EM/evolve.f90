	module M_EVOLVE
	implicit none
	
	contains
	
	
	subroutine evolve(unew,uold,dx,dt,time,x,y,z,Nv,N)
	use m_derivs
	use m_boundary
	use m_rhs
	use m_pars, only : STMETH, disip, iEx,iEy,iEz,iff,iBx,iBy,iBz,igg
	implicit none
	
	real*8, dimension(:,:,:,:), intent(inout):: unew
	real*8, dimension(:,:,:,:), intent(in):: uold
	real*8, dimension(:,:,:):: x,y,z	
	real*8 dx, dt, time
	integer :: Nv,N
	
!keep in some variables the rhs
	real*8, dimension(:,:,:,:), allocatable :: dxU,dyU,dzU,urk1,urk2,urk3,urk4,disu
	real*8 EPS, factor2,factor3
	integer i,j
	logical :: ltrace

!for debug
	ltrace = .false.
	
!allocate memmory needed to take all steps
!dxU and dyU will keep the derivatives of the field
!urk1,2,3,4 are to keep the intermidiate rhs values for the RK update
! the derivatives, rhs, boundaries have their on files where the routines
! are kept.

	allocate(dxU(Nv,N,N,N),dyU(Nv,N,N,N),dzU(Nv,N,N,N), &
        &        urk1(Nv,N,N,N),urk2(Nv,N,N,N),&
	&        urk3(Nv,N,N,N),urk4(Nv,N,N,N),&
        &        disu(Nv,N,N,N))
	
! We'll do the rk 3rd or 4th  order . The procedure is basically
! the same 3 times, (i) evaluate derivs; (ii) evaluate rhs;
! (iii) get the new value
	
	EPS = disip 

	if(STMETH.eq.3) then
	factor2 = 0.75
	else if(STMETH.eq.4) then
	factor2 = 0.5
	else
	print*, 'nocoded!'
	end if	
!!!!!!!!!!!!!!!!!!! FIRST STEP	
		
	do i = 1,Nv
	  call derivs(uold(i,:,:,:),dxU(i,:,:,:),dx,N,1)
	  call derivs(uold(i,:,:,:),dyU(i,:,:,:),dx,N,2)
	  call derivs(uold(i,:,:,:),dzU(i,:,:,:),dx,N,3)
	  call dissip(uold(i,:,:,:),disu(i,:,:,:),dx,N)
	end do

	call rhs(uold,dxU,dyU,dzU,urk1,dx,N,x)

	
        if(ltrace) print*, 'call boun'
	call boundary(uold,dxU,dyU,dzU,urk1,dx,time,x)
	
	urk1 = urk1 + EPS*disu
	
	unew = uold + dt * 0.5 * urk1
	

!!!!!!!!!!!!!!!!!!! SECOND STEP	

!calculate derivs	
	do i = 1,Nv
	  call derivs(unew(i,:,:,:),dxU(i,:,:,:),dx,N,1)
	  call derivs(unew(i,:,:,:),dyU(i,:,:,:),dx,N,2)
	  call derivs(unew(i,:,:,:),dzU(i,:,:,:),dx,N,3)
	  call dissip(unew(i,:,:,:),disu(i,:,:,:),dx,N)	  
	end do
	
	call rhs(unew,dxU,dyU,dzU,urk2,dx,N,x)
	
	call boundary(unew,dxU,dyU,dzU, urk2,dx,time,x)
	
	urk2 = urk2 + EPS*disu

	unew = uold + dt * factor2 * urk2	
	
!!!!!!!!!!!!!!!!!!! THIRD STEP	

!calculate derivs	
	do i = 1,Nv
	  call derivs(unew(i,:,:,:),dxU(i,:,:,:),dx,N,1)
	  call derivs(unew(i,:,:,:),dyU(i,:,:,:),dx,N,2)
	  call derivs(unew(i,:,:,:),dzU(i,:,:,:),dx,N,3)
	  call dissip(unew(i,:,:,:),disu(i,:,:,:),dx,N)		  
	end do
	
	
	call rhs(unew,dxU,dyU,dzU,urk3,dx,N,x)
	
	call boundary(unew,dxU,dyU,dzU, urk3,dx,time,x)
	
	urk3 = urk3 + EPS*disu

	unew = uold + dt * urk3

!!!!!!!!!!!!!!!!!!! FOURTH STEP	

   if (STMETH.eq.4) THEN	
!calculate derivs	
	do i = 1,Nv
	  call derivs(unew(i,:,:,:),dxU(i,:,:,:),dx,N,1)
	  call derivs(unew(i,:,:,:),dyU(i,:,:,:),dx,N,2)
	  call derivs(unew(i,:,:,:),dzU(i,:,:,:),dx,N,3)
	  call dissip(unew(i,:,:,:),disu(i,:,:,:),dx,N)	  
	end do
	
	call rhs(unew,dxU,dyU,dzU,urk4,dx,N,x)
	
	call boundary(unew,dxU,dyU,dzU, urk4,dx,time,x)
	
	urk4 = urk4 + EPS*disu

        unew = uold + dt/6. * (urk1+2.*urk2+2.*urk3+urk4)
   else 	
	unew = uold + dt/9. * (2.*urk1+3.*urk2+4.*urk3)		
   end if


!deallocate memmory
	deallocate(dxU,dyU,dzU,urk1,urk2,urk3,urk4,disu)

	
	end subroutine evolve
	
	end module M_EVOLVE
