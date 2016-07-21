	MODULE M_BOUNDARY
	use m_pars, only:  sl, sr, derorder, c, s,a
	
	implicit none
	
	contains
	
	
	subroutine boundary(u,dxu,dyu,dzU,uboun,dx,time,x)
	use m_pars, only:  n, sl, sr, derorder, c, s,a, gamma
	implicit none	
	real*8, dimension(:,:,:,:), intent(in) :: u,dxu,dyu,dzU
	real*8, dimension(:,:,:,:), intent(inout) :: uboun
	real*8, dimension(:,:,:), intent(inout) :: x
	real*8 dx, time, factor
	real*8 Vp, Vm, vel,nn, cl2
	integer i,j


	real*8 :: UU(7,7), TT(7,7),cx(4),cy(4),GNx(4),GNy(4)
	real*8 :: dp, dm, B2, Bn, cs, s1, Bt, Rp, BM, Rm, BV
	real*8 :: vn,pn,vt,Bti,dl, sigma
	real*8 :: vbn,pbn,detl,F1,F2,F3,F4,dets,CF,DF,EF,FF
	real*8 :: tempbase1, tempbase2, NNt, Ly,vy,vx,Bx,By
	real*8 :: rho, p, phi, pb, tempcheck, FCP(2), AA(7,7)
	real*8 :: nx, ny, RHS(7), lam(7), templam(7), tmap, temp
	integer :: ilam, jlam, li, lj ,lk,lm
        real*8 :: Rproj(7)

!set periodic boundaries

! 	uboun(:,:,1) = 0.5*(uboun(:,:,1)+uboun(:,:,N))
! 	uboun(:,:,N) = uboun(:,:,1)	
! 	
! 	uboun(:,1,:) = 0.5*(uboun(:,1,:)+uboun(:,N,:))
! 	uboun(:,N,:) = uboun(:,1,:)


	end subroutine boundary
	
	
	end module M_BOUNDARY
	
	
	
