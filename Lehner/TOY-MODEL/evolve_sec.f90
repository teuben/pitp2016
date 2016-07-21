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
        real*8  :: dx, dt, sigma, timeb, timec, aa, alpha, beta, eta
        real*8  :: onehalf, onethird, twothird, onesixth, onefourth
        real*8  :: oneeigth ,onetuelveth
        real*8  :: at21,at31,at32,at41,at42,at43,omt1,omt2,omt3,omt4
        real*8  :: a11,a21,a22,a31,a32,a33,a41,a42,a43,a44,om1,om2,om3,om4

	
!keep in some variables the rhs
	real*8, dimension(:,:),allocatable:: rk1,rk2,rk3,rk4
	real*8, dimension(:,:),allocatable:: ru1,ru2,ru3,ru4
	real*8, dimension(:,:),allocatable:: du,diss,u1,u2,u3,u4
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
	allocate(rk1(dim,Nx),rk2(dim,Nx),rk3(dim,Nx),rk4(dim,Nx))
	allocate(u1(dim,Nx),u2(dim,Nx),u3(dim,Nx),u4(dim,Nx))
	allocate(ru1(dim,Nx),ru2(dim,Nx),ru3(dim,Nx),ru4(dim,Nx))

       !  remember u(1)=V and u(2)=W

        !-----------------------------------------------------------
        ! DIFFERENT TYPES OF IMEX RK--------------------------------
        !-----------------------------------------------------------
 
        !  IMEX-SSP2 (3,3,2) stiffly accurate scheme 
        if (runge_kutta .EQ. 1) then
          at21 = 0.0d0
          at31 = 0.0d0
          at32 = onehalf
          at41 = 0.0d0
          at42 = onehalf
          at43 = onehalf
          omt1 = 0.0d0
          omt2 = onethird
          omt3 = onethird
          omt4 = onethird

          a11 = 0.0d0
          a21 = 0.0d0
          a22 = onefourth
          a31 = 0.0d0
          a32 = 0.0d0
          a33 = onefourth
          a41 = 0.0d0
          a42 = onethird
          a43 = onethird
          a44 = onethird
          om1 = 0.0d0
          om2 = onethird
          om3 = onethird
          om4 = onethird

        !  IMEX-SSP3 (3,3,2) L-stable scheme
        else if (runge_kutta .EQ. 2) then
          at21 = 0.0d0
          at31 = 0.0d0
          at32 = 1.0d0
          at41 = 0.0d0
          at42 = onefourth
          at43 = onefourth
          omt1 = 0.0d0
          omt2 = onesixth
          omt3 = onesixth
          omt4 = twothird

          alpha = 1.0d0 - 1.0d0/dsqrt(2.0d0)
          a11 = 0.0d0
          a21 = 0.0d0
          a22 = alpha
          a31 = 0.0d0
          a32 = 1.0d0 - 2.0d0*alpha
          a33 = alpha
          a41 = 0.0d0
          a42 = onehalf - alpha
          a43 = 0.0d0
          a44 = alpha
          om1 = 0.0d0
          om2 = onesixth
          om3 = onesixth
          om4 = twothird

        !  IMEX-SSP3 (4,3,3) L-stable scheme
        else if (runge_kutta .EQ. 3) then
          at21 = 0.0d0
          at31 = 0.0d0
          at32 = 1.0d0
          at41 = 0.0d0
          at42 = onefourth
          at43 = onefourth
          omt1 = 0.0d0
          omt2 = onesixth
          omt3 = onesixth
          omt4 = twothird

          alpha = 0.24169426078821
          beta  = 0.06042356519705
          eta   = 0.12915286960590

          a11 = alpha
          a21 =-alpha
          a22 = alpha
          a31 = 0.0d0
          a32 = 1.0d0 - alpha
          a33 = alpha
          a41 = beta
          a42 = eta
          a43 = onehalf - beta - eta - alpha
          a44 = alpha
          om1 = 0.0d0
          om2 = onesixth
          om3 = onesixth
          om4 = twothird

        !  IMEX-SSP3 (4,3,3) L-stable scheme & almost Boscarini 
        else if (runge_kutta .EQ. 4) then
          at21 = 0.0d0
          at31 = 0.0d0
          at32 = 1.0d0
          at41 = 0.0d0
          at42 = onefourth
          at43 = onefourth
          omt1 = 0.0d0
          omt2 = onesixth
          omt3 = onesixth
          omt4 = twothird

          a11 = onehalf
          a21 =-onehalf
          a22 = onehalf
          a31 = onesixth
          a32 = onehalf
          a33 = onethird
          a41 = onetuelveth
          a42 = onetuelveth
          a43 = 0.0d0
          a44 = onethird
          om1 = 0.0d0
          om2 = onesixth
          om3 = onesixth
          om4 = twothird

        !  IMEX-SSP3 (4,3,3) L-stable scheme simpler
        else if (runge_kutta .EQ. 5) then
          at21 = 0.0d0
          at31 = 0.0d0
          at32 = 1.0d0
          at41 = 0.0d0
          at42 = onefourth
          at43 = onefourth
          omt1 = 0.0d0
          omt2 = onesixth
          omt3 = onesixth
          omt4 = twothird

          a11 = onehalf
          a21 =-onehalf
          a22 = onehalf
          a31 = 0.0d0
          a32 = onehalf
          a33 = onehalf
          a41 = oneeigth
          a42 = 0.0d0
          a43 = oneeigth
          a44 = onehalf
          om1 = 0.0d0
          om2 = onesixth
          om3 = onesixth
          om4 = twothird

        !  IMEX-SSP3 (4,4,2) L-stable scheme stiffly accurate
        else if (runge_kutta .EQ. 6) then
          at21 = 1.0d0
          at31 = 2.0d0/9.0d0
          at32 = 1.0d0/9.0d0
          at41 = 0.0d0
          at42 = 0.0d0
          at43 = 1.0d0
          omt1 = 0.0d0
          omt2 = 0.0d0
          omt3 = 3.0d0/4.0d0
          omt4 = 1.0d0/4.0d0

          a11 = 1.0d0/4.0d0
          a21 = 3.0d0/4.0d0
          a22 = 1.0d0/4.0d0
          a31 = 4.0d0/27.0d0
          a32 = -7.0d0/108.d0
          a33 = 1.0d0/4.0d0
          a41 = 0.0d0
          a42 = 0.0d0
          a43 = 3.0d0/4.0d0
          a44 = 1.0d0/4.0d0
          om1 = 0.0d0
          om2 = 0.0d0
          om3 = 3.0d0/4.0d0
          om4 = 1.0d0/4.0d0

        !  IMEX-SSP3 (4,4,2) L-stable scheme stiffly accurate II
        else if (runge_kutta .EQ. 7) then
          at21 = 1.0d0/3.0d0
          at31 = 0.0d0
          at32 = 1.0d0/3.0d0
          at41 = 0.0d0
          at42 = 0.0d0
          at43 = 1.0d0
          omt1 = 0.0d0
          omt2 = 0.0d0
          omt3 = 3.0d0/4.0d0
          omt4 = 1.0d0/4.0d0

          a11 = 1.0d0/2.0d0
          a21 = 1.0d0/12.0d0
          a22 = 1.0d0/4.0d0
          a31 = -13.0d0/24.0d0
          a32 = 3.0d0/8.d0
          a33 = 1.0d0/2.0d0
          a41 = 0.0d0
          a42 = 1.0d0/4.0d0
          a43 = 1.0d0/2.0d0
          a44 = 1.0d0/4.0d0
          om1 = 0.0d0
          om2 = 1.0d0/4.0d0
          om3 = 1.0d0/2.0d0
          om4 = 1.0d0/4.0d0

        else
          print*, "This RK is not implemented, runge_kutta=",runge_kutta
          stop
        end if
	
!	Fren4  = 1.0d0 / ( 1.0d0 + onefourth * cond * dt)
!	Fren3  = 1.0d0 / ( 1.0d0 + onethird * cond * dt)

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


        ! FINAL STEP-----------------------------	
	!-----------------------------------------

        do i = pos,nx
          unew(:,i) = uold(:,i) &
          & + dt*( omt1*rk1(:,i) + omt2*rk2(:,i) + omt3*rk3(:,i) + omt4*rk4(:,i)) &
          & + dt*( om1*ru1(:,i)  +  om2*ru2(:,i) +  om3*ru3(:,i) +  om4*ru4(:,i) )
        end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!deallocate memmory
	deallocate(du,rk1,rk2,rk3,rk4,diss)
	deallocate(u1,u2,u3,u4)
	deallocate(ru1,ru2,ru3,ru4)
	
	end subroutine evolve
	
	
	end module M_EVOLVE
