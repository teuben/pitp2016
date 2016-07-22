	program main
		
	use m_initial
	use m_evolve
	use m_pars 
	
	implicit none

	
!define variables, of an undetermined size yet

	real*8, dimension(:,:), allocatable :: unew, uold
	real*8, dimension(:), allocatable :: x, mout, conduc
	
!variables for the rest
	real*8 :: dx, dt, Rmax, maskf, pp
	integer :: i,j, gft_out_full,gft_out_brief, ret, eb, ebo
        character(32) :: formato
        formato = '(F10.4, E16.8, E16.8)'
	

!read pars
	call readpars
   	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
!allocate memmory
	allocate(uold(dim,Nx),unew(dim,nx),x(nx),mout(nx), conduc(nx))
!define grid
	dx = (Rout)/dble(Nx-1)
	dt = cfl * dx
	do i=1, nx
	  x(i) = dble(i-1)*dx
!	  conduc(i) = 1.0 * cond * exp( -((x(i) - 0.3)/0.08)**2  )
	  conduc(i) = 1.0 * cond 
	end do

        mout = 1
        eb = 1
	
!intialize
	time = 0.0
	call initial(uold,x,Nx)
        print*, time, maxval(abs(uold(1,:))),maxval(abs(uold(2,:)))
!        ret = gft_out_full('conduc',time, nx-1, 'x', 1, x,conduc(1:nx-1))
!        ret = gft_out_full('V',time, nx-1, 'x', 1, x,uold(1,1:nx-1)) 
!        ret = gft_out_full('W',time, nx-1, 'x', 1, x,uold(2,1:nx-1)) 

!save for analysis
        open(unit = 100, file="Us.dat")
!        write(100,*) "#time=",time
!        do j=1,Nx
!          write(100,formato), x(j), uold(1,j), uold(2,j)
!        end do
			
!time loop	
	do i = 1, Nt

          call evolve(unew,uold,conduc,dx,dt,eb,x)

!update fileds and time	
          uold = unew
	  time = time + dt

	  ! compute the new conductivity
!          conduc = 1.0 * cond * exp( -((x - 0.3)/0.08)**2  )
	
!every now and then do some output
	  if(mod(i,freq).eq.0) then
	    !will call here for output
	    print*, time, maxval(abs(uold(1,:))),maxval(abs(uold(2,:)))
 !           ret = gft_out_full('conduc',time, nx-1, 'x', 1, x,conduc(1:nx-1))
!            ret = gft_out_full('V',time, nx-1, 'x', 1, x,uold(1,1:nx-1))
!            ret = gft_out_full('W',time, nx-1, 'x', 1, x,uold(2,1:nx-1))

!            write(100,*) 
!            write(100,*) "#time=",time
!            do j=1,Nx
!              write(100,formato), x(j), uold(1,j), uold(2,j)
!            end do
	  end if	
	
	end do

        do j=1,Nx
          write(100,formato), x(j), uold(1,j), uold(2,j)
        end do
        close(100)	

	end program
