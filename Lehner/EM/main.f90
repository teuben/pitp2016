!program to integrate source free maxwell's equations
!very simple, uniform 3d grid
!also, artificially setting periodic boundary conditions

	program main
	use M_initial
	use M_evolve
	use M_pars
	use M_derivs
	
	implicit none
	
	INTERFACE
	subroutine norm3d(norm,infield,N)
	implicit none
	real*8, dimension(:,:,:):: infield
	real*8 :: norm
	integer N

	end subroutine norm3d

	END INTERFACE

!U will be conservative vars
!P will be primitive vars
	real*8, dimension(:,:,:,:), allocatable :: Unew, Uold	
	real*8, dimension(:,:,:), allocatable :: x, y, z,dxBx,dyBy,dzBz
	real*8 :: dx,dy,dz, dt, time, divBnorm
	integer :: ITE, i,j, Nv, Nsteps
        integer :: gft_out_full,gft_out_brief, ret, coun, ic, INFO
	logical :: ltrace	


!set for debug or not
	ltrace = .false.

!read pars
        call readpars


!set the max number of eqns for now
	Nv = 8


!assume same extent for both coordinates
	allocate(Unew(Nv,N,N,N), Uold(Nv,N,N,N), &
     &           x(N,N,N),y(N,N,N),z(N,N,N), &
     &           dxBx(N,N,N),dyBy(N,N,N),dzBz(N,N,N))
	
		
!define coords and related	
	
	dx= 10./(N-1)
	dy = dx
	dz = dx
	dt = cfl * dx
	
	do i = 1, N
	x(i,:,:) = -5.+ (i-1)*dx
	y(:,i,:) = -5.+ (i-1)*dy
	z(:,:,i) = -5.+ (i-1)*dz
	end do
	
!define initial data
	call initial(Uold,dx,dy,dz,x,y,z,Nv,N)
        time = 0.0d0	
	Unew = Uold

	if(ltrace) print*, 'about to integrate'
	
	   do ITE=1, Nt
	    !integrate

        if(ltrace) print*, 'call evolve'
	    call evolve(unew,uold,dx,dt,time,x,y,z,Nv,N)	
        if(ltrace) print*, 'back from evolve'

          if (mod(ITE-1,freq).eq.0) then
          print*, 'call out', time
!calculate vorticity
	  call derivs(uold(iBx,:,:,:),dxBx,dx,N,1)
	  call derivs(uold(iBy,:,:,:),dyBy,dx,N,2)
	  call derivs(uold(iBz,:,:,:),dzBz,dx,N,3)

!          ret = gft_out_brief('constraint',time, (/n,n,n/), 3, dxBx+dyBy+dzBz)  
!          ret = gft_out_brief('Bx',time, (/n,n,n/), 3, uold(iBx,:,:,:))
!          ret = gft_out_brief('By',time, (/n,n,n/), 3, uold(iBy,:,:,:))
!          ret = gft_out_brief('Bz',time, (/n,n,n/), 3, uold(iBz,:,:,:))

	  call norm3d(divBnorm,dxBx+dyBy+dzBz,N)

	  write(101,*), time, divBnorm
 
        if(ltrace) print*, 'end call out'
	  end if

	   !update the dime and save the new values
	   !in the old array. 

            time = time + dt
            Uold = Unew

	   end do


	end

	subroutine norm3d(norm,infield,N)
	implicit none
	real*8, dimension(:,:,:):: infield
	real*8 :: norm
	integer N
	integer :: i,j,k

	norm = 0.0d0

	do i =1, N
	 do j=1, N
	  do k=1, N
	   norm = norm + infield(i,j,k)**2
	  end do
	 end do
	end do

	norm = sqrt(norm/(N*N*N))
	
	end subroutine norm3d
