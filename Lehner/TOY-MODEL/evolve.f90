	module M_EVOLVE
	implicit none
	
	contains
	
	
	subroutine evolve(unew,uold,conduc,dx,dt, eb,x)
	use m_derivs
	use m_boundary
	use m_rhs
	use m_pars, only : disip, dim, c, runge_kutta
	implicit none
	
	integer eb, method
	real*8, dimension(:,:), intent(inout):: unew
	real*8, dimension(:,:), intent(in):: uold
	real*8, dimension(:) :: x, conduc
        real*8  :: dx, dt, sigma, timeb, timec, aa, alpha, delta, eta
        real*8  :: onehalf, onethird, twothird, onesixth, onefourth
        real*8  :: oneeigth ,onetuelveth
        real*8  :: at21,at31,at32,at41,at42,at43,at51,at52,at53,at54
        real*8  :: a11,a21,a22,a31,a32,a33,a41,a42,a43,a44,a51,a52,a53,a54,a55
        real*8  :: omt1,omt2,omt3,omt4,omt5,om1,om2,om3,om4,om5
	
!keep in some variables the rhs
	real*8, dimension(:,:),allocatable:: rk1,rk2,rk3,rk4,rk5
	real*8, dimension(:,:),allocatable:: ru1,ru2,ru3,ru4,ru5
	real*8, dimension(:,:),allocatable:: du,diss,u1,u2,u3,u4,u5
	integer i,j

        pos = 1 
        onehalf  = 1.0d0/2.0d0
        onethird = 1.0d0/3.0d0
        twothird = 2.0d0/3.0d0
        onesixth = 1.0d0/6.0d0
        oneeigth = 1.0d0/8.0d0
        onetuelveth = 1.0d0/12.0d0
        onefourth = 1.0d0/4.0d0
	sigma = -disip
	
	
!allocate memmory
	allocate(du(dim,Nx),diss(dim,Nx))
	allocate(rk1(dim,Nx),rk2(dim,Nx),rk3(dim,Nx),rk4(dim,Nx),rk5(dim,Nx))
	allocate(u1(dim,Nx),u2(dim,Nx),u3(dim,Nx),u4(dim,Nx),u5(dim,Nx))
	allocate(ru1(dim,Nx),ru2(dim,Nx),ru3(dim,Nx),ru4(dim,Nx),ru5(dim,Nx))

       !  remember u(1)=V and u(2)=W

        !-----------------------------------------------------------
        ! DIFFERENT TYPES OF IMEX RK--------------------------------
        !-----------------------------------------------------------
 
        !  IMEX-SSP3 (5,3,2) L-stable scheme alpha=2/3
        if (runge_kutta .EQ. 1) then

          at21 = 0.0d0
          at31 = 0.0d0
          at32 = 1.0d0
          at41 = 0.0d0
          at42 = 1.0d0/4.0d0
          at43 = 1.0d0/4.0d0
          at51 = 0.0d0
          at52 = 1.0d0/6.0d0
          at53 = 1.0d0/6.0d0
          at54 = 2.0d0/3.0d0
          omt1 = 0.0d0
          omt2 = 1.0d0/6.0d0
          omt3 = 1.0d0/6.0d0
          omt4 = 2.0d0/3.0d0
          omt5 = 0.0d0

          alpha = 4.0d0/6.0d0

          a11 = alpha
          a21 =-alpha
          a22 = alpha
          a31 = 0.0d0
          a32 = 1.0d0 - alpha
          a33 = alpha
          a41 = (2.0d0*alpha*alpha - 1.0d0 + 2.0d0*alpha)/(8.0d0*alpha)
          a42 =-(4.0d0*alpha*alpha - 1.0d0)/(8.0d0*alpha)
          a43 = -3.0d0*alpha/4.0d0 + 1.0d0/4.0d0
          a44 = alpha
          a51 = 0.0d0
          a52 = 1.0d0/6.0d0
          a53 = 0.0d0
          a54 = 2.0d0/3.0d0
          a55 = 1.0d0/6.0d0
          om1 = 0.0d0
          om2 = 1.0d0/6.0d0
          om3 = 0.0d0
          om4 = 2.0d0/3.0d0
          om5 = 1.0d0/6.0d0

        !  IMEX-SSP3 (5,3,2) L-stable scheme alpha=5/6
        else if (runge_kutta .EQ. 2) then

          at21 = 0.0d0
          at31 = 0.0d0
          at32 = 1.0d0
          at41 = 0.0d0
          at42 = 1.0d0/4.0d0
          at43 = 1.0d0/4.0d0
          at51 = 0.0d0
          at52 = 1.0d0/6.0d0
          at53 = 1.0d0/6.0d0
          at54 = 2.0d0/3.0d0
          omt1 = 0.0d0
          omt2 = 1.0d0/6.0d0
          omt3 = 1.0d0/6.0d0
          omt4 = 2.0d0/3.0d0
          omt5 = 0.0d0

          alpha = 5.0d0/6.0d0

          a11 = alpha
          a21 =-alpha
          a22 = alpha
          a31 = 0.0d0
          a32 = 1.0d0 - alpha
          a33 = alpha
          a41 = (2.0d0*alpha*alpha - 1.0d0 + 2.0d0*alpha)/(8.0d0*alpha)
          a42 =-(4.0d0*alpha*alpha - 1.0d0)/(8.0d0*alpha)
          a43 = -3.0d0*alpha/4.0d0 + 1.0d0/4.0d0
          a44 = alpha
          a51 = 0.0d0
          a52 = 1.0d0/6.0d0
          a53 = 0.0d0
          a54 = 2.0d0/3.0d0
          a55 = 1.0d0/6.0d0
          om1 = 0.0d0
          om2 = 1.0d0/6.0d0
          om3 = 0.0d0
          om4 = 2.0d0/3.0d0
          om5 = 1.0d0/6.0d0

        !  IMEX-SSP3 (5,3,2) L-stable scheme alpha=5/6
        else if (runge_kutta .EQ. 3) then

          at21 = 0.0d0
          at31 = 0.0d0
          at32 = 1.0d0
          at41 = 0.0d0
          at42 = 1.0d0/4.0d0
          at43 = 1.0d0/4.0d0
          at51 = 0.0d0
          at52 = 1.0d0/6.0d0
          at53 = 1.0d0/6.0d0
          at54 = 2.0d0/3.0d0
          omt1 = 0.0d0
          omt2 = 1.0d0/6.0d0
          omt3 = 1.0d0/6.0d0
          omt4 = 2.0d0/3.0d0
          omt5 = 0.0d0

          a11 = 0.5d0
          a21 =-0.5d0
          a22 = 0.5d0
          a31 = 0.0d0
          a32 = 0.5d0
          a33 = 0.5d0
          a41 = 1.0d0/8.0d0
          a42 = 0.0d0
          a43 =-1.0d0/8.0d0
          a44 = 0.5d0
          a51 = 0.0d0
          a52 = 1.0d0/6.0d0
          a53 = 1.0d0/6.0d0
          a54 = 2.0d0/3.0d0
          a55 = 0.0d0
          om1 = 0.0d0
          om2 = 1.0d0/6.0d0
          om3 = 1.0d0/6.0d0
          om4 = 2.0d0/3.0d0
          om5 = 0.0d0

        else
          print*, "This RK is not implemented, runge_kutta=",runge_kutta
          stop
        end if
	
        ! FIRST STEP-----------------------------	
	!-----------------------------------------
	! compute the explicit part----------
        do i=pos,nx	
	  u1(:,i) =  uold(:,i)
	end do 

	! compute the implicit part----------
        do i=pos,nx	
   	  aa   = a11 * conduc(i) * dt
     	  u1(1,i) =( u1(1,i) + aa * c * u1(2,i) ) / (1.0d0 + aa)
	end do
	
        ! compute the rhs and R------
	do i = 1, dim
	  call derivs(u1(i,:),du(i,:),dx,Nx,eb)
          call dissip(u1(i,:),diss(i,:),dx,Nx,eb)
	end do
	
	call rhs(u1,du,rk1,ru1,conduc,dx,Nx,eb,x)
	rk1 = rk1 + sigma*diss

        call boundary(u1,du,rk1,dx)

        ! SECOND STEP-----------------------------
	!-----------------------------------------
	! compute W* and V*----------
        do i=pos,nx	
	  u2(:,i) = uold(:,i) + at21 * dt * rk1(:,i) + a21  * dt * ru1(:,i)
	end do

	! compute the implicit part----------
        do i=pos,nx	
   	  aa   = a22 * conduc(i) * dt
     	  u2(1,i) =( u2(1,i) + aa * c * u2(2,i) ) / (1.0d0 + aa)
	end do

        ! compute the rhs and R------
	do i = 1, dim
	  call derivs(u2(i,:),du(i,:),dx,Nx,eb)
          call dissip(u2(i,:),diss(i,:),dx,Nx,eb)
	end do
	
	call rhs(u2,du,rk2,ru2,conduc,dx,Nx,eb,x)
	rk2 = rk2 + sigma*diss	

        call boundary(u2,du,rk2,dx)

        ! THIRD STEP-----------------------------
	!-----------------------------------------
	! compute W* and V*----------
        do i=pos,nx	
          u3(:,i) = uold(:,i) + dt * ( at31 * rk1(:,i) + at32 * rk2(:,i) ) &
                            & + dt * (  a31 * ru1(:,i) +  a32 * ru2(:,i) )
	end do    

	! compute the implicit part----------
        do i=pos,nx	
   	  aa   = a33 * conduc(i) * dt
     	  u3(1,i) =( u3(1,i) + aa * c * u3(2,i) ) / (1.0d0 + aa)
	end do

        ! compute the rhs and R------
	do i = 1, dim
	  call derivs(u3(i,:),du(i,:),dx,Nx,eb)
          call dissip(u3(i,:),diss(i,:),dx,Nx,eb)
	end do
	
	call rhs(u3,du,rk3,ru3,conduc,dx,Nx,eb,x)
	rk3 = rk3 + sigma*diss

        call boundary(u3,du,rk3,dx)

        ! FOURTH STEP-----------------------------
	!-----------------------------------------
	! compute W* and V*----------
        do i=pos,nx	
          u4(:,i) = uold(:,i) &
                  & + dt * ( at41 * rk1(:,i) + at42 * rk2(:,i) + at43 * rk3(:,i)) &
                  & + dt * (  a41 * ru1(:,i) +  a42 * ru2(:,i) +  a43 * ru3(:,i))
	end do    

	! compute the implicit part----------
        do i=pos,nx	
   	  aa   = a44 * conduc(i) * dt
     	  u4(1,i) =( u4(1,i) + aa * c * u4(2,i) ) / (1.0d0 + aa)
	end do

        ! compute the rhs and R------
	do i = 1, dim
	  call derivs(u4(i,:),du(i,:),dx,Nx,eb)
          call dissip(u4(i,:),diss(i,:),dx,Nx,eb)
	end do
	
	call rhs(u4,du,rk4,ru4,conduc,dx,Nx,eb,x)
	rk4 = rk4 + sigma*diss

        call boundary(u4,du,rk4,dx)

        ! FIFTH STEP------------------------------
	!-----------------------------------------
	! compute W* and V*----------
        do i=pos,nx	
          u5(:,i) = uold(:,i) &
            & + dt * ( at51 * rk1(:,i) + at52 * rk2(:,i) + at53 * rk3(:,i) + at54 * rk4(:,i)) &
            & + dt * (  a51 * ru1(:,i) +  a52 * ru2(:,i) +  a53 * ru3(:,i) +  a54 * ru4(:,i))
	end do    

	! compute the implicit part----------
        do i=pos,nx	
   	  aa   = a55 * conduc(i) * dt
     	  u5(1,i) =( u5(1,i) + aa * c * u5(2,i) ) / (1.0d0 + aa)
	end do


        ! FINAL STEP-----------------------------	
	!-----------------------------------------

        unew = u5

!        do i = pos,nx
!          unew(:,i) = uold(:,i) &
!          & + dt*( omt1*rk1(:,i) + omt2*rk2(:,i) + omt3*rk3(:,i) + omt4*rk4(:,i)) &
!          & + dt*( om1*ru1(:,i)  +  om2*ru2(:,i) +  om3*ru3(:,i) +  om4*ru4(:,i) )
!        end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!deallocate memmory
	deallocate(du,rk1,rk2,rk3,rk4,rk5,diss)
	deallocate(u1,u2,u3,u4,u5)
	deallocate(ru1,ru2,ru3,ru4,ru5)
	
	end subroutine evolve
	
	
	end module M_EVOLVE
