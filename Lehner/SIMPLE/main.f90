	program main

!load modules of things needed
		
	use m_initial
	use m_evolve
	use m_pars 
	
	implicit none

	
!define variables, of an undetermined size yet
! using 9, too many for this example but it is a
! stripped down version
! x radial grid

	real*8, dimension(:,:), allocatable :: unew, uold
	real*8, dimension(:), allocatable :: x
	
!variables for the rest
	real*8 :: dx, dt, Rmax, time
	integer :: i,j, gft_out_full,gft_out_brief, ret
	
	integer :: m1, m2, m3

!read pars
	call readpars
   	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
!allocate memmory
	allocate(uold(9,Nx),unew(9,nx),x(nx))	
!define grid
	dx = (-Rmin+Rout)/(Nx-1)
	dt = cfl * dx
	do i=1, nx
	x(i) = Rmin + (i-1)*dx
	end do

	m1 = int((40-Rmin)/dx) + 1
        m2 = int((70-Rmin)/dx) + 1
        m3 = int((100-Rmin)/dx) + 1

!intialize
	time = 0.0
	call initial(uold,x,Nx)
	
		
!time loop	
	do i = 1, Nt
	
	call evolve(unew,uold,dx,dt,x)
	
!update fileds and time	
	uold = unew   	
	time = time + dt
	
!every now and then do some output
	  if(mod(i,freq).eq.0) then
	  !will call here for output
	  print*, 'outputting' 
!  	   ret = gft_out_full('phi',time, nx, 'x', 1, x, uold(1,:))
!          ret = gft_out_full('g',time, nx, 'x', 1, x, uold(2,:))
!          ret = gft_out_full('f',time, nx, 'x', 1, x, uold(3,:))

!for output, ascii ----
	do j=1,Nx
	write(101,*) x(j), uold(1,j)
	write(102,*) x(j), uold(2,j)
	write(103,*) x(j), uold(3,j)
	end do
	write(101,*) 
	write(102,*) 
	write(103,*) 


	write(101,*) 
	write(102,*) 
	write(103,*) 

	  end if

	
	end do
	
	end
