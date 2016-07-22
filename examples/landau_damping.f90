program landau_damping

  ! program solves dg/dt + v * d/dx (g + phi)
  ! phi = int dv exp(-v**2) * g

  implicit none

  !-------------- input parameters ---------------------!
  
  ! Te/Ti
  real, parameter :: tau = 1.0
  
  ! xgrid parameters
  ! nx is total number of x grid points
  integer, parameter :: nx = 20
  real, parameter :: Lx = 1.0

  ! vgrid parameters
  ! nvp is number of positive parallel velocities sampled
  integer, parameter :: nvp = 100
  real, parameter :: Lv = 3.0

  ! time grid paramters
  integer, parameter :: nstep = 1000
  real, parameter :: dt = 0.001
  integer, parameter :: nwrite = 10

  ! k-space parameters
  ! k = 2*pi*kint
!  integer :: kint = 1
  
  !------------------------------------------------------!
  
  ! grid spacing in x and v
  real :: dx, dv

  ! total number of parallel velocities sampled
  integer :: nv

  ! keeps track of time variable
  real :: time
  
  ! the constant pi
  real :: pi
  
  ! g(x,v)
  real, dimension (:,:), allocatable :: gxv
  ! phi(x)
  real, dimension (:), allocatable :: phi
  ! grid in x
  real, dimension (:), allocatable :: xgrid
  ! grid in v
  real, dimension (:), allocatable :: vgrid

  ! gk(m)
!  complex, dimension (:), allocatable :: gkm
  
  integer :: phi_unit = 101, gxv_unit = 102, phi2_unit = 103
  
  ! define pi for later use
  pi = 2.*acos(0.)
!  k = 2.*pi*kint
  
  call init_io
  
  ! allocates and populates xgrid and computes dx
  call get_xgrid
  ! allocates and populates vgrid and computes dv  
  call get_vgrid

  ! allocates and initializes g(x,v)
  call init_gxv
  ! allocates and initializes phi(x)
  call init_phi
  ! explicit in time, grid-based in x and v solver
  call simple_solve

  ! allocates and initializes g(k,m)
!  call init_gkm
  ! implicit in time, spectral in x and v solver
!  call advanced_solve
  
  call finish_io

  deallocate (xgrid, vgrid)
  deallocate (gxv)
  deallocate (phi)
  
contains

  subroutine init_io

    implicit none

    open (phi_unit, file='landau_damping.phi', status='replace')
    open (phi2_unit, file='landau_damping.phi2', status='replace')
    open (gxv_unit, file='landau_damping.gxv', status='replace')
    
  end subroutine init_io

  subroutine finish_io

    implicit none

    close (phi_unit)
    close (phi2_unit)
    close (gxv_unit)
    
  end subroutine finish_io
  
  subroutine get_xgrid

    implicit none

    integer :: i
    
    if (.not.allocated(xgrid)) allocate (xgrid(nx))

    dx = Lx/(nx-1)
    do i = 1, nx
       xgrid(i) = (i-1)*dx
    end do
    
  end subroutine get_xgrid

  subroutine get_vgrid
    
    implicit none

    integer :: i
    
    if (.not.allocated(vgrid)) allocate (vgrid(-nvp:nvp))

    nv = 2*nvp+1
    dv = Lv/nvp
    do i = -nvp, nvp
       vgrid(i) = i*dv
    end do
    
  end subroutine get_vgrid

  subroutine init_gxv

    implicit none

    if (.not.allocated(gxv)) allocate (gxv(nx,-nvp:nvp))

    gxv = spread(cos(2.*pi*xgrid),2,nv)
    
  end subroutine init_gxv

  ! subroutine init_gkm

  !   implicit none

  !   if (.not.allocated(gkm)) allocate (gkm(0:nm-1))

  !   gkm(0) = 0.5
  !   gkm(1:) = 0.0
    
  ! end subroutine init_gkm
  
  subroutine init_phi

    implicit none

    if (.not.allocated(phi)) allocate (phi(nx))

    call get_phi_grid
    
  end subroutine init_phi

  subroutine get_phi_grid

    implicit none

    integer :: i
    
    do i = 1, nx
       call integrate_v (gxv(i,:),phi(i))
    end do
    phi = tau*phi
    
  end subroutine get_phi_grid
  
  subroutine integrate_v (f, tot)

    implicit none

    real, dimension (-nvp:), intent (in) :: f
    real, intent (out) :: tot

    integer :: i

    tot = 0.
    do i = -nvp, nvp
       tot = tot + dv*exp(-vgrid(i)**2)*f(i)/sqrt(pi)
    end do
    
  end subroutine integrate_v
  
  ! simple_solve solves Landau damping problem
  ! with explicit time-stepping (forward Euler)
  ! upwind differencing in x
  ! simple grid-based integration in v
  subroutine simple_solve

    implicit none

    integer :: istep

    time = 0.0
    call write_simple (0)
    
    do istep = 1, nstep
       call update_gxv
       call get_phi_grid
       time = time + dt
       call write_simple (istep)
    end do
    
  end subroutine simple_solve

  subroutine update_gxv

    implicit none

    integer :: i
    
    ! special case of v = 0 -- dg/dt = 0
    ! so no need to do anything

    ! v < 0
    do i = 1, nx-1
       gxv(i,:-1) = gxv(i,:-1) + dt*vgrid(:-1)/dx &
            * (gxv(i,:-1) - gxv(i+1,:-1) &
            + phi(i) - phi(i+1))
    end do
    ! need a boundary condition on x
    ! use periodicity
    gxv(nx,:-1) = gxv(1,:-1)
    
    ! v > 0
    do i = nx, 2, -1
       gxv(i,1:) = gxv(i,1:) + dt*vgrid(1:)/dx &
            * (gxv(i-1,1:) - gxv(i,1:) &
            + phi(i-1) - phi(i))
    end do
    ! need a boundary condition on x
    ! use periodicity
    gxv(1,1:) = gxv(nx,1:)
    
  end subroutine update_gxv

  subroutine write_simple (istep)

    implicit none

    integer, intent (in) :: istep
    
    integer :: i, j

    if (mod(istep,nwrite)==0) then
       do j = -nvp, nvp
          do i = 1, nx
             write (gxv_unit,*) time, xgrid(i), vgrid(j), gxv(i,j)
          end do
          write (gxv_unit,*)
       end do
       write (gxv_unit,*)
       write (gxv_unit,*)
       
       do i = 1, nx
          write (phi_unit,*) time, xgrid(i), phi(i)
       end do
       write (phi_unit,*)
    end if
       
    write (phi2_unit,*) time, sum(phi**2)

    write (*,*) 'time = ', time, '|phi|**2 = ', sum(phi**2)
    
  end subroutine write_simple

!  subroutine advanced_solve

!    implicit none

!    integer :: istep

!    time = 0.0
!    call write_advanced(0)
    
!    do istep = 1, nstep
!       call update_gkm
!       time = time + dt
!       call write_advanced (istep)
!    end do
    
!  end subroutine advanced_solve

  subroutine update_gkm

    implicit none
    
  end subroutine update_gkm

  subroutine write_advanced

    implicit none
    
  end subroutine write_advanced
  
  ! solves system Ax = b for x (which is returned as sol)
  ! inputs are aa, bb, and cc (the elements to the left, center, and right
  ! of diagonal in tridiagonal matrix A)
  ! and sol=b as the rhs of the linear system Ax = b
  subroutine tridag (aa, bb, cc, sol)

    implicit none

    real, dimension (:), intent (in) :: aa, bb, cc
    real, dimension (:), intent (in out) :: sol

    integer :: ix, npts
    real :: bet

    real, dimension (:), allocatable :: gam
    
    npts = size(aa)
    allocate (gam(npts))

    bet = bb(1)
    sol(1) = sol(1)/bet

    do ix = 2, npts
       gam(ix) = cc(ix-1)/bet
       bet = bb(ix) - aa(ix)*gam(ix)
       if (bet == 0.0) write (*,*) 'tridiagonal solve failed'
       sol(ix) = (sol(ix)-aa(ix)*sol(ix-1))/bet
    end do

    do ix = npts-1, 1, -1
       sol(ix) = sol(ix) - gam(ix+1)*sol(ix+1)
    end do

    deallocate (gam)

  end subroutine tridag
  
end program landau_damping
