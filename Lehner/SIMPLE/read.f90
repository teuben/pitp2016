        program main
        implicit none
        real*8, dimension(:), allocatable :: gxx, x
        integer :: i,j, gft_read_full,gft_out_full, ret, &
	&	   gft_read_shape, nx
	character*32 :: cname
	real*8 :: time
	
!first make sure we know of the size of the grid	
	ret=gft_read_shape('grr',1,nx)
	
!allocate enough space
	allocate(x(nx),gxx(nx))		
		
	ret=gft_read_full('grr',1,nx,cname,1,time,x,gxx)
	
!write to check we get the same
        ret = gft_out_full('grr2',time, nx, 'x', 1, x, gxx)

	end 
