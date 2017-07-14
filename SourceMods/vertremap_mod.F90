!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Begin GPU remap module  !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vertremap_mod

  !**************************************************************************************
  !
  !  Purpose:
  !        Construct sub-grid-scale polynomials using piecewise spline method with
  !        monotone filters.
  !
  !  References: PCM - Zerroukat et al., Q.J.R. Meteorol. Soc., 2005. (ZWS2005QJR)
  !              PSM - Zerroukat et al., Int. J. Numer. Meth. Fluids, 2005. (ZWS2005IJMF)
  !
  !**************************************************************************************

  use shr_kind_mod,           only: r8=>shr_kind_r8
  use dimensions_mod,         only: np,nlev,qsize,nlevp,npsq,nc
  use hybvcoord_mod,          only: hvcoord_t
  use element_mod,            only: element_t
  use fvm_control_volume_mod, only: fvm_struct
  use perf_mod,               only: t_startf, t_stopf ! _EXTERNAL
  use parallel_mod,           only: parallel_t
  use cam_abortutils,         only: endrun
  !use control_mod,            only: vert_remap_q_alg

  public remap1                  ! remap any field, splines, monotone
  public remap1_nofilter         ! remap any field, splines, no filter
! todo: tweak interface to match remap1 above, rename remap1_ppm:
  public remap_q_ppm             ! remap state%Q, PPM, monotone

  contains

!=======================================================================================================!

subroutine remap1(Qdp,nx,qstart,qstop,qsize,dp1,dp2,hybrid)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass, monotone on Q=Qdp/dp
  !
  use hybrid_mod, only : hybrid_t, get_loop_ranges, config_thread_region
  use thread_mod, only : tracer_num_threads
  implicit none
  integer, intent(in) :: nx,qstart,qstop,qsize
  real (kind=r8), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=r8), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  type (hybrid_t), optional :: hybrid
  ! ========================
  ! Local Variables
  ! ========================

  type (hybrid_t) :: hybridnew
  real (kind=r8), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real (kind=r8), dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp
  real (kind=r8)  :: f_xm,level1,level2,level4,level5, &
       peaks_min,peaks_max,tmp_cal,xm,xm_d,zv1,zv2, &
       zero = 0._r8,one = 1._r8,tiny = 1.e-12_r8,qmax = 1.e50_r8
  integer :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, &
       lt1,lt2,lt3,t1,t2,t3,t4,tm,tp,i,ilev,j,jk,k,q
  integer :: qbeg, qend,vert_remap_q_alg=1
  logical :: abort=.false.

  if (vert_remap_q_alg == 1 .or. vert_remap_q_alg == 2) then
    call t_startf('remap_Q_ppm')
    if ( present(hybrid) ) then
      !$OMP PARALLEL NUM_THREADS(tracer_num_threads), DEFAULT(SHARED), PRIVATE(hybridnew,qbeg,qend)
      hybridnew = config_thread_region(hybrid,'tracer')
      call get_loop_ranges(hybridnew, qbeg=qbeg, qend=qend)
      call remap_Q_ppm(qdp,nx,qbeg,qend,qsize,dp1,dp2,1)
      !$OMP END PARALLEL
    else
      call remap_Q_ppm(qdp,nx,qstart,qstop,qsize,dp1,dp2,1)
    endif
    call t_stopf('remap_Q_ppm')
    return
  endif

  !  write(6,*) 'YOU ARE DOINGS SPLINES REMAP WHICH IS NOT STRICTLY MONOTONE! ABORT'
  !  abort = .true.                                                                 
  do q=qstart,qstop
    do i=1,nx
      do j=1,nx

        z1c(1)=0 ! source grid
        z2c(1)=0 ! target grid
        do k=1,nlev
          z1c(k+1)=z1c(k)+dp1(i,j,k)
          z2c(k+1)=z2c(k)+dp2(i,j,k)
        enddo

        zv(1)=0
        do k=1,nlev
          Qcol(k)=Qdp(i,j,k,q)!  *(z1c(k+1)-z1c(k)) input is mass
          zv(k+1) = zv(k)+Qcol(k)
        enddo

        if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001_r8) then
          write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
          write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
          write(6,*) 'DATA FOR MODEL LEVELS'
          write(6,*) 'PLEVMODEL=',z2c(nlev+1)
          write(6,*) 'PLEV     =',z1c(nlev+1)
          write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
          abort=.true.
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! quadratic splies with UK met office monotonicity constraints  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        zkr  = 99
        ilev = 2
        zkr(1) = 1
        zkr(nlev+1) = nlev
        kloop: do k = 2,nlev
          do jk = ilev,nlev+1
            if (z1c(jk).ge.z2c(k)) then
              ilev      = jk
              zkr(k)   = jk-1
              cycle kloop
            endif
          enddo
        enddo kloop

        zgam  = (z2c(1:nlev+1)-z1c(zkr)) / (z1c(zkr+1)-z1c(zkr))
        zgam(1)      = 0.0_r8
        zgam(nlev+1) = 1.0_r8
        zhdp = z1c(2:nlev+1)-z1c(1:nlev)


        h = 1/zhdp
        zarg = Qcol * h
        rhs = 0
        lower_diag = 0
        diag = 0
        upper_diag = 0

        rhs(1)=3*zarg(1)
        rhs(2:nlev) = 3*(zarg(2:nlev)*h(2:nlev) + zarg(1:nlev-1)*h(1:nlev-1))
        rhs(nlev+1)=3*zarg(nlev)

        lower_diag(1)=1
        lower_diag(2:nlev) = h(1:nlev-1)
        lower_diag(nlev+1)=1

        diag(1)=2
        diag(2:nlev) = 2*(h(2:nlev) + h(1:nlev-1))
        diag(nlev+1)=2

        upper_diag(1)=1
        upper_diag(2:nlev) = h(2:nlev)
        upper_diag(nlev+1)=0

        q_diag(1)=-upper_diag(1)/diag(1)
        rhs(1)= rhs(1)/diag(1)

        do k=2,nlev+1
          tmp_cal    =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
          q_diag(k) = -upper_diag(k)*tmp_cal
          rhs(k) =  (rhs(k)-lower_diag(k)*rhs(k-1))*tmp_cal
        enddo
        do k=nlev,1,-1
          rhs(k)=rhs(k)+q_diag(k)*rhs(k+1)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!  monotonicity modifications  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        filter_code = 0
        dy(1:nlev-1) = zarg(2:nlev)-zarg(1:nlev-1)
        dy(nlev) = dy(nlev-1)

        dy = merge(zero, dy, abs(dy) < tiny )

        do k=1,nlev
          im1=MAX(1,k-1)
          im2=MAX(1,k-2)
          im3=MAX(1,k-3)
          ip1=MIN(nlev,k+1)
          t1 = merge(1,0,(zarg(k)-rhs(k))*(rhs(k)-zarg(im1)) >= 0)
          t2 = merge(1,0,dy(im2)*(rhs(k)-zarg(im1)) > 0 .AND. dy(im2)*dy(im3) > 0 &
               .AND. dy(k)*dy(ip1) > 0 .AND. dy(im2)*dy(k) < 0 )
          t3 = merge(1,0,ABS(rhs(k)-zarg(im1)) > ABS(rhs(k)-zarg(k)))

          filter_code(k) = merge(0,1,t1+t2 > 0)
          rhs(k) = (1-filter_code(k))*rhs(k)+filter_code(k)*(t3*zarg(k)+(1-t3)*zarg(im1))
          filter_code(im1) = MAX(filter_code(im1),filter_code(k))
        enddo

        rhs = merge(qmax,rhs,rhs > qmax)
        rhs = merge(zero,rhs,rhs < zero)

        za0 = rhs(1:nlev)
        za1 = -4*rhs(1:nlev) - 2*rhs(2:nlev+1) + 6*zarg
        za2 =  3*rhs(1:nlev) + 3*rhs(2:nlev+1) - 6*zarg

        dy(1:nlev) = rhs(2:nlev+1)-rhs(1:nlev)
        dy = merge(zero, dy, abs(dy) < tiny )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Compute the 3 quadratic spline coeffients {za0, za1, za2}                 !!
        !! knowing the quadratic spline parameters {rho_left,rho_right,zarg}         !!
        !! Zerroukat et.al., Q.J.R. Meteorol. Soc., Vol. 128, pp. 2801-2820 (2002).   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        h = rhs(2:nlev+1)

        do k=1,nlev
          xm_d = merge(one,2*za2(k),abs(za2(k)) < tiny)
          xm = merge(zero,-za1(k)/xm_d, abs(za2(k)) < tiny)
          f_xm = za0(k) + za1(k)*xm + za2(k)*xm**2

          t1 = merge(1,0,ABS(za2(k)) > tiny)
          t2 = merge(1,0,xm <= zero .OR. xm >= 1)
          t3 = merge(1,0,za2(k) > zero)
          t4 = merge(1,0,za2(k) < zero)
          tm = merge(1,0,t1*((1-t2)+t3) .EQ. 2)
          tp = merge(1,0,t1*((1-t2)+(1-t3)+t4) .EQ. 3)

          peaks=0
          peaks = merge(-1,peaks,tm .EQ. 1)
          peaks = merge(+1,peaks,tp .EQ. 1)
          peaks_min = merge(f_xm,MIN(za0(k),za0(k)+za1(k)+za2(k)),tm .EQ. 1)
          peaks_max = merge(f_xm,MAX(za0(k),za0(k)+za1(k)+za2(k)),tp .EQ. 1)

          im1=MAX(1,k-1)
          im2=MAX(1,k-2)
          ip1=MIN(nlev,k+1)
          ip2=MIN(nlev,k+2)

          t1 = merge(abs(peaks),0,(dy(im2)*dy(im1) <= tiny) .OR. &
               (dy(ip1)*dy(ip2) <= tiny) .OR. (dy(im1)*dy(ip1) >= tiny) .OR. &
               (dy(im1)*float(peaks) <= tiny))

          filter_code(k) = merge(1,t1+(1-t1)*filter_code(k),(rhs(k) >= qmax) .OR. &
               (rhs(k) <= zero) .OR. (peaks_max > qmax) .OR. (peaks_min < tiny))

          if (filter_code(k) > 0) then
            level1 = rhs(k)
            level2 = (2*rhs(k)+h(k))/3
            level4 = (1/3_r8)*rhs(k)+2*(1/3_r8)*h(k)
            level5 = h(k)

            t1 = merge(1,0,h(k) >= rhs(k))
            t2 = merge(1,0,zarg(k) <= level1 .OR.  zarg(k) >= level5)
            t3 = merge(1,0,zarg(k) >  level1 .AND. zarg(k) <  level2)
            t4 = merge(1,0,zarg(k) >  level4 .AND. zarg(k) <  level5)

            lt1 = t1*t2
            lt2 = t1*(1-t2+t3)
            lt3 = t1*(1-t2+1-t3+t4)

            za0(k) = merge(zarg(k),za0(k),lt1 .EQ. 1)
            za1(k) = merge(zero,za1(k),lt1 .EQ. 1)
            za2(k) = merge(zero,za2(k),lt1 .EQ. 1)

            za0(k) = merge(rhs(k),za0(k),lt2 .EQ. 2)
            za1(k) = merge(zero,za1(k),lt2 .EQ. 2)
            za2(k) = merge(3*(zarg(k)-rhs(k)),za2(k),lt2 .EQ. 2)

            za0(k) = merge(-2*h(k)+3*zarg(k),za0(k),lt3 .EQ. 3)
            za1(k) = merge(+6*h(k)-6*zarg(k),za1(k),lt3 .EQ. 3)
            za2(k) = merge(-3*h(k)+3*zarg(k),za2(k),lt3 .EQ. 3)

            t2 = merge(1,0,zarg(k) >= level1 .OR.  zarg(k) <= level5)
            t3 = merge(1,0,zarg(k) <  level1 .AND. zarg(k) >  level2)
            t4 = merge(1,0,zarg(k) <  level4 .AND. zarg(k) >  level5)

            lt1 = (1-t1)*t2
            lt2 = (1-t1)*(1-t2+t3)
            lt3 = (1-t1)*(1-t2+1-t3+t4)

            za0(k) = merge(zarg(k),za0(k),lt1 .EQ. 1)
            za1(k) = merge(zero,za1(k),lt1 .EQ. 1)
            za2(k) = merge(zero,za2(k),lt1 .EQ. 1)

            za0(k) = merge(rhs(k),za0(k),lt2 .EQ. 2)
            za1(k) = merge(zero,za1(k),lt2 .EQ. 2)
            za2(k) = merge(3*(zarg(k)-rhs(k)),za2(k),lt2 .EQ. 2)

            za0(k) = merge(-2*h(k)+3*zarg(k),za0(k),lt3 .EQ. 3)
            za1(k) = merge(+6*h(k)-6*zarg(k),za1(k),lt3 .EQ. 3)
            za2(k) = merge(-3*h(k)+3*zarg(k),za2(k),lt3 .EQ. 3)
          endif
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! start iteration from top to bottom of atmosphere !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        zv1 = 0
        do k=1,nlev
          if (zgam(k+1)>1_r8) then
            WRITE(*,*) 'r not in [0:1]', zgam(k+1)
            abort=.true.
          endif
          zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
               (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
          Qdp(i,j,k,q) = (zv2 - zv1) ! / (z2c(k+1)-z2c(k) ) dont convert back to mixing ratio
          zv1 = zv2
        enddo
      enddo
    enddo
  enddo ! q loop
  if (abort) then
    call endrun('Bad levels in remap1.  usually CFL violatioin')
  end if

end subroutine remap1

subroutine remap1_nofilter(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  implicit none
  integer, intent(in) :: nx,qsize
  real (kind=r8), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=r8), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! ========================
  ! Local Variables
  ! ========================

  real (kind=r8), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real (kind=r8), dimension(nlev)      :: h,Qcol,za0,za1,za2,zarg,zhdp
  real (kind=r8)  :: tmp_cal,zv1,zv2
  integer :: zkr(nlev+1),i,ilev,j,jk,k,q
  logical :: abort=.false.
  !   call t_startf('remap1_nofilter')

#if (defined COLUMN_OPENMP)
  !$omp parallel do num_threads(tracer_num_threads) &
  !$omp    private(q,i,j,z1c,z2c,zv,k,Qcol,zkr,ilev) &
  !$omp    private(jk,zgam,zhdp,h,zarg,rhs,lower_diag,diag,upper_diag,q_diag,tmp_cal) &
  !$omp    private(za0,za1,za2) &
  !$omp    private(ip2,zv1,zv2)
#endif
  do q=1,qsize
    do i=1,nx
      do j=1,nx

        z1c(1)=0 ! source grid
        z2c(1)=0 ! target grid
        do k=1,nlev
          z1c(k+1)=z1c(k)+dp1(i,j,k)
          z2c(k+1)=z2c(k)+dp2(i,j,k)
        enddo

        zv(1)=0
        do k=1,nlev
          Qcol(k)=Qdp(i,j,k,q)!  *(z1c(k+1)-z1c(k)) input is mass
          zv(k+1) = zv(k)+Qcol(k)
        enddo

        if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001_r8) then
          write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
          write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
          write(6,*) 'DATA FOR MODEL LEVELS'
          write(6,*) 'PLEVMODEL=',z2c(nlev+1)
          write(6,*) 'PLEV     =',z1c(nlev+1)
          write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
          abort=.true.
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! quadratic splies with UK met office monotonicity constraints  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        zkr  = 99
        ilev = 2
        zkr(1) = 1
        zkr(nlev+1) = nlev
        kloop: do k = 2,nlev
          do jk = ilev,nlev+1
            if (z1c(jk).ge.z2c(k)) then
              ilev      = jk
              zkr(k)   = jk-1
              cycle kloop
            endif
          enddo
        enddo kloop

        zgam  = (z2c(1:nlev+1)-z1c(zkr)) / (z1c(zkr+1)-z1c(zkr))
        zgam(1)      = 0.0_r8
        zgam(nlev+1) = 1.0_r8
        zhdp = z1c(2:nlev+1)-z1c(1:nlev)


        h = 1/zhdp
        zarg = Qcol * h
        rhs = 0
        lower_diag = 0
        diag = 0
        upper_diag = 0

        rhs(1)=3*zarg(1)
        rhs(2:nlev) = 3*(zarg(2:nlev)*h(2:nlev) + zarg(1:nlev-1)*h(1:nlev-1))
        rhs(nlev+1)=3*zarg(nlev)

        lower_diag(1)=1
        lower_diag(2:nlev) = h(1:nlev-1)
        lower_diag(nlev+1)=1

        diag(1)=2
        diag(2:nlev) = 2*(h(2:nlev) + h(1:nlev-1))
        diag(nlev+1)=2

        upper_diag(1)=1
        upper_diag(2:nlev) = h(2:nlev)
        upper_diag(nlev+1)=0

        q_diag(1)=-upper_diag(1)/diag(1)
        rhs(1)= rhs(1)/diag(1)

        do k=2,nlev+1
          tmp_cal    =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
          q_diag(k) = -upper_diag(k)*tmp_cal
          rhs(k) =  (rhs(k)-lower_diag(k)*rhs(k-1))*tmp_cal
        enddo
        do k=nlev,1,-1
          rhs(k)=rhs(k)+q_diag(k)*rhs(k+1)
        enddo

        za0 = rhs(1:nlev)
        za1 = -4*rhs(1:nlev) - 2*rhs(2:nlev+1) + 6*zarg
        za2 =  3*rhs(1:nlev) + 3*rhs(2:nlev+1) - 6*zarg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! start iteration from top to bottom of atmosphere !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        zv1 = 0
        do k=1,nlev
          if (zgam(k+1)>1_r8) then
            WRITE(*,*) 'r not in [0:1]', zgam(k+1)
            abort=.true.
          endif
          zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
               (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
          Qdp(i,j,k,q) = (zv2 - zv1) ! / (z2c(k+1)-z2c(k) ) dont convert back to mixing ratio
          zv1 = zv2
        enddo
      enddo
    enddo
  enddo ! q loop
  if (abort) then
    call endrun('Bad levels in remap1_nofilter.  usually CFL violatioin')
  end if
end subroutine remap1_nofilter

!=============================================================================!


!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
!using PPM instead of splines.
subroutine remap_Q_ppm(Qdp,nx,qstart,qstop,qsize,dp1,dp2,field)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  !use control_mod, only        : vert_remap_q_alg
  use spmd_utils            , only: masterproc
  implicit none
  integer,intent(in) :: nx,qstart,qstop,qsize
  real (kind=r8), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=r8), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  integer, intent(in) :: field !1 = ttmp remap, 2 = velocities
  ! Local Variables
  integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
  real(kind=r8), dimension(  1-gs:nlev+2 ) :: pio    !Pressure at interfaces for old grid
  real(kind=r8), dimension(  1-gs:nlev+2 ) :: pio_c  !Pressure at cell for old grid 
  real(kind=r8), dimension( 4  ) :: x_i,f_i  !Pressure at cell for old grid 
  real(kind=r8), dimension(  1-gs:nlev+2 ) :: lag    !Lagrangian interpolator 
  real(kind=r8), dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
  real(kind=r8), dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
  real(kind=r8), dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
  real(kind=r8), dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
  real(kind=r8), dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
  real(kind=r8), dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
  real(kind=r8), dimension(       nlev   ) :: z1, z2
  real(kind=r8) :: ppmdx(10,0:nlev+1)  !grid spacings
  real(kind=r8) :: massn1, massn2,pio_l
  integer :: i, j, k, q, kk,kid(nlev),vert_remap_q_alg=1,m,unitn=8,x_p,extrap,ao_extrap
  character(len=256):: filename


  do j = 1 , nx
    do i = 1 , nx

      pin(1)=0
      pio(1)=0
      do k=1,nlev
         dpn(k)=dp2(i,j,k)
         dpo(k)=dp1(i,j,k)
         pin(k+1)=pin(k)+dpn(k)
         pio(k+1)=pio(k)+dpo(k)
      enddo

      pio(nlev+2) = pio(nlev+1) + 1._r8  !This is here to allow an entire block of k threads to run in the remapping phase.
                                      !It makes sure there's an old interface value below the domain that is larger.
      pin(nlev+1) = pio(nlev+1)       !The total mass in a column does not change.
                                      !Therefore, the pressure of that mass cannot either.
      !Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence.
      do k = 1 , gs
         dpo(1   -k) = dpo(       k)  
         dpo(nlev+k) = dpo(nlev+1-k)
      enddo

      extrap = 1 ! 0 = mirrored, 1 = linear, 2 = lagrangian
      !defining the pressure inside cell to extrapolate ao later

      if(extrap==1)then
      !LINEAR:
      !extrapolating pio
             pio(0) = pio(2)+ (pio(1)-pio(2)) + dpo(1)
             pio(-1) = pio(1) + ( 2. )*(pio(0)-pio(1)) + dpo(0)
 
            ! pio(nlev+1) = pio(nlev-1) + ( (nlev+1)-(nlev-1))/( (nlev)-(nlev-1) )*(pio(nlev)-pio(nlev-1))+dpo(nlev)
            ! pio(nlev+2) = pio(nlev) + ( (nlev+2)-(nlev))/( (nlev+1)-nlev)*(pio(nlev+1)-pio(nlev))+dpo(nlev+1)

      !defining the pressure inside cell to extrapolate ao later
             do k=0,nlev+1
                pio_c(k) = (pio(k+1)+pio(k))/2.0_r8
             end do

             pio_c(-1) = pio_c(1) + ( 2. )*(pio_c(0)-pio_c(1))
             !pio_c(nlev+2) = pio_c(nlev) + ( ((nlev+2) - (nlev))/( (nlev+1) - (nlev))  )*(pio_c(nlev+1)-pio_c(nlev)) 
        
      !LAGRANGIAN:
      elseif(extrap==2)then

      x_i(1:4) =  (/2, 4, 6, 11/) !taking 4 indexes to interpolate
      f_i(1:4) =  (/pio(2), pio(4), pio(6), pio(11)/) !and it's value iny-axis

      do x_p = 0,-1,-1
       pio_l = 0.0d0
       do m = 1,4
         lag = 1.0_r8
         do kk = 1,4
                if(kk/=m)then
                    lag(m) = lag(m)*( x_p  - x_i(kk))/( x_i(m) - x_i(kk) )
                end if
         enddo
         pio_l = pio_l + f_i(m)*lag(m)
       enddo
      pio(x_p) = pio_l
      enddo




      x_i(1:4) =  (/2, 4, 6, 11/) !taking 4 indexes to interpolate
      f_i(1:4) =  (/pio_c(2), pio_c(4), pio_c(6), pio_c(11)/) !and it's value in y-axis

      do x_p = 0,-1,-1
       pio_l = 0.0d0
       do m = 1,4
         lag = 1.0_r8
         do kk = 1,4
                if(kk/=m)then
                    lag(m) = lag(m)*( x_p  - x_i(kk))/( x_i(m) - x_i(kk) )
                end if
         enddo
         pio_l = pio_l + f_i(m)*lag(m)
       enddo
      pio_c(x_p) = pio_l
      enddo
     end if


     if(masterproc .and. field==1 )then
           write (filename, '("pio_extrap.dat")' )
           open(unitn, file=trim(filename),status='replace' )
           do k=-1,nlev+2
               write(unitn,*) pio(k)
           enddo
           close(unitn)
     endif

      !Compute remapping intervals once for all tracers. Find the old grid cell index in which the
      !k-th new cell interface resides. Then integrate from the bottom of that old cell to the new
      !interface location. In practice, the grid never deforms past one cell, so the search can be
      !simplified by this. Also, the interval of integration is usually of magnitude close to zero
      !or close to dpo because of minimial deformation.
      !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so
      !I set them equal to each other.
      do k = 1 , nlev
        kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
        !Find the index of the old grid cell in which this new cell's bottom interface resides.
        do while ( pio(kk) <= pin(k+1) )
          kk = kk + 1
        enddo
        kk = kk - 1                   !kk is now the cell index we're integrating over.
        if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds.
                                      !Top bounds match anyway, so doesn't matter what coefficients are used
        kid(k) = kk                   !Save for reuse
        z1(k) = -0.5_R8                !This remapping assumes we're starting from the left interface of an old grid cell
                                      !In fact, we're usually integrating very little or almost all of the cell in question
        z2(k) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5_r8 ) / dpo(kk)  !PPM interpolants are normalized to an independent
                                                                        !coordinate domain [-0.5,0.5].
      enddo

      !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the
      !limiting. So anything that depends only on the grid is pre-computed outside the tracer loop.
      ppmdx(:,:) = compute_ppm_grids( dpo )

      !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
      !mass accumulation
      do q = qstart, qstop
        !Accumulate the old mass up to old grid cell interface locations to simplify integration
        !during remapping. Also, divide out the grid spacing so we're working with actual tracer
        !values and can conserve mass. The option for ifndef ZEROHORZ I believe is there to ensure
        !tracer consistency for an initially uniform field. I copied it from the old remap routine.
        masso(1) = 0._r8

        do k = 1 , nlev
          ao(k) = Qdp(i,j,k,q)
          masso(k+1) = masso(k) + ao(k) !Accumulate the old mass. This will simplify the remapping
          ao(k) = ao(k) / dpo(k)        !Divide out the old grid spacing because we want the tracer mixing ratio, not mass.
        enddo

        extrap = 1
        ao_extrap = 1

        !Fill in ghost values. Ignored if vert_remap_q_alg == 2
        do k = 1 , gs
          ao(1   -k) = ao(       k)
          ao(nlev+k) = ao(nlev+1-k)
        enddo

        !LINEAR:
        if(extrap==1)then
         !extrapolating ao (scalar quant.):
         !left:
!         if(pio_c(1) - pio_c(2) /= 0) then
!                        ao(0) = ao(2) + (( pio_c(0) - pio_c(2) )/( pio_c(1) - pio_c(2)  ))*( ao(1 )- ao(2) )  
                        ao(0) = ao(2) + (( 0 - 2 )/( 1 - 2  ))*( ao(1)- ao(2) )  
!         else
!                ao(0) = ao(1)
!         endif

!         if( pio_c(0) - pio_c(1) /=0)then
!                        ao(-1) = ao(1) + (( pio_c(-1) - pio_c(1) )/( pio_c(0) - pio_c(1) ))*( ao(0)- ao(1) )   
                        ao(-1) = ao(1) + (( -1 - 1 )/( 0 - 1 ))*( ao(0)- ao(1) )   
!         else
!                 ao(-1) = ao(0)
!         endif

         !if(pio(nlev) - pio(nlev-1)/=0)then
         !right:
!                 ao(nlev+1) = ao(nlev-1) + ( pio_c(nlev+1) - pio_c(nlev-1) )/(pio_c(nlev) - pio_c(nlev-1) )*&
!                       ( ao(nlev) - ao(nlev-1) )
         !else
         !       ao(nlev+1) = ao(nlev)
         !endif

!         if(pio_c(nlev+1) - pio_c(nlev)/=0)then
!                 ao(nlev+2) = ao(nlev) + ( pio_c(nlev+2) - pio_c(nlev))/(pio_c(nlev+1) - pio_c(nlev) )*&
!                       ( ao(nlev+1) - ao(nlev) )
!         else
!                 ao(nlev+2) = ao(nlev+1)
!         endif
       endif

       if(ao_extrap==1)then
       !LAGRANGIAN:
      
       x_i(1:4) =  (/1., 4., 6., 15./)                 !taking 4 indexes to interpolate ! 2-4-6-11
!       x_i(1:4) =  (/pio_c(2), pio_c(4), pio_c(6), pio_c(11)/)                 !taking 4 indexes to interpolate 
       f_i(1:4) =  (/ao(1), ao(4), ao(6), ao(15)/) !and it's value in y-axis

       
      do x_p = 0,-1,-1
       pio_l = 0.0d0    
       do m = 1,4
         lag = 1.0_r8
         do kk = 1,4
                if(m/=kk)then
                    lag(m) = lag(m)*( x_p  - x_i(kk))/( x_i(m) - x_i(kk) ) !x_i indices
!                    lag(m) = lag(m)*( pio_c(x_p)  - x_i(kk))/( x_i(m) - x_i(kk) )
                end if
         enddo
         pio_l = pio_l + f_i(m)*lag(m)
       enddo      
       ao(x_p) = pio_l
      enddo


       x_i(1:4) =  (/nlev, nlev-4, nlev-6, nlev-11/)                
!       x_i(1:4) =  (/pio_c(nlev), pio_c(nlev-4), pio_c(nlev-6), pio_c(nlev-11)/) 
       f_i(1:4) =  (/ao(nlev), ao(nlev-4), ao(nlev-6), ao(nlev-11)/) 

      do x_p = nlev+1,nlev+2
       pio_l = 0.0d0    
       do m = 1,4
         lag = 1.0_r8
         do kk = 1,4
                if(m/=kk)then
                    lag(m) = lag(m)*( x_p  - x_i(kk))/( x_i(m) - x_i(kk) ) !x_i indices
!                    lag(m) = lag(m)*( pio_c(x_p)  - x_i(kk))/( x_i(m) - x_i(kk))
                end if
         enddo
         pio_l = pio_l + f_i(m)*lag(m)
       enddo      
!       ao(x_p) = pio_l
      enddo
      end if


      if(masterproc .and. field==1 )then
            write (filename, '("ao_extrap.dat")' )
            open(unitn, file=trim(filename),status='replace') !'old',position='append' )
            do k=-1,nlev+2
                if(k>0 .and. k<nlev+1)then
                write(unitn,*) ao(k),Qdp(1,1,k,1)
                else
                write(unitn,*) ao(k)
                endif
            enddo
            close(unitn)
       endif


        !Compute monotonic and conservative PPM reconstruction over every cell
        coefs(:,:) = compute_ppm( ao , ppmdx )
        !Compute tracer values on the new grid by integrating from the old cell bottom to the new
        !cell interface to form a new grid mass accumulation. Taking the difference between
        !accumulation at successive interfaces gives the mass inside each cell. Since Qdp is
        !supposed to hold the full mass this needs no normalization.
        massn1 = 0._r8
        do k = 1 , nlev
          kk = kid(k)
          massn2 = masso(kk) + integrate_parabola( coefs(:,kk) , z1(k) , z2(k) ) * dpo(kk)
          Qdp(i,j,k,q) = massn2 - massn1
          massn1 = massn2
        enddo
      enddo
    enddo
  enddo
! call t_stopf('remap_Q_ppm')
end subroutine remap_Q_ppm

subroutine remap2(Qdp,nx,qstart,qstop,qsize,dp_1,dp_2,rmp_kind,hybrid)

  !Wrapper for the vertical remapping techniques explored by M. Valera during
  !the SiParCS internship, Summer 2017.
  !Layout based on previous code (remap1())/

  !============
  ! rmp_kind(s):
  ! 1.- filtered (remap1_filter())
  ! 2.- ppm (remap_Q_ppm)
  ! 3.- mom (remapping_core_h() using method from vertical_remap() line 1131
  ! -> initialize_remapping())
  ! 4.- no filter (remap1_nofilter())
  !============
 
  ! remap field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass, monotone on Q=Qdp/dp
  !
  use hybrid_mod, only : hybrid_t, get_loop_ranges, config_thread_region
  use thread_mod, only : tracer_num_threads
  use MOM_remapping

  implicit none
  integer, intent(in) :: nx,qstart,qstop,qsize
  real (kind=r8), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=r8), intent(in) :: dp_1(nx,nx,nlev),dp_2(nx,nx,nlev)
  type (hybrid_t), optional :: hybrid
  integer, intent(in) :: rmp_kind
  ! ========================
  ! Local Variables
  ! ========================

  type (hybrid_t) :: hybridnew
  real (kind=r8), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real (kind=r8), dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp
  real (kind=r8)  :: f_xm,level1,level2,level4,level5, &
       peaks_min,peaks_max,tmp_cal,xm,xm_d,zv1,zv2, &
       zero = 0._r8,one = 1._r8,tiny = 1.e-12_r8,qmax = 1.e50_r8
  integer :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, &
       lt1,lt2,lt3,t1,t2,t3,t4,tm,tp,i,ilev,j,jk,k,q
  integer :: qbeg, qend,vert_remap_q_alg=1
  logical :: abort=.false.

  real (kind=r8), dimension(np,np,nlev)  :: dp_1_inv!,dp_2_inv
  real, dimension(nlev) :: dp_1_mom,qdp_1_mom,dp_2_mom,qdp_2_mom

  dp_1_inv = 1._r8/dp_1
!  dp_2_inv = 1._r8/dp_2

  if(rmp_kind==1)then
            call remap1(Qdp,np,1,1,1,dp_1,dp_2) !T_rmp*dp_moist
  elseif(rmp_kind==2)then
            call remap_Q_ppm(Qdp,np,1,1,1,dp_1,dp_2,1)
  elseif(rmp_kind==3)then

          do i=1,np
                do j=1,np

                    dp_1_mom = dp_1(i,j,:)
                    dp_2_mom = dp_2(i,j,:)
                    qdp_1_mom = Qdp(i,j,:,qsize)*dp_1_inv(i,j,:) !qsize = np1 ??
                    qdp_2_mom = 0._r8

                    call remapping_core_h(CSP,nlev,dp_1_mom,qdp_1_mom,nlev,dp_2_mom,qdp_2_mom)

                    Qdp(i,j,:,qsize) = qdp_2_mom
                enddo
          enddo

          Qdp(:,:,:,qsize) = Qdp(:,:,:,qsize)*dp_2   

  else
            call remap1_nofilter(Qdp,np,2,dp_1,dp_2)
  end if

end subroutine remap2

!=======================================================================================================!


!THis compute grid-based coefficients from Collela & Woodward 1984.
function compute_ppm_grids( dx )   result(rslt)
  !use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=r8), intent(in) :: dx(-1:nlev+2)  !grid spacings
  real(kind=r8)             :: rslt(10,0:nlev+1)  !grid spacings
  integer :: j
  integer :: indB, indE,vert_remap_q_alg=1

  !Calculate grid-based coefficients for stage 1 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    rslt( 1,j) = dx(j) / ( dx(j-1) + dx(j) + dx(j+1) )
    rslt( 2,j) = ( 2._r8*dx(j-1) + dx(j) ) / ( dx(j+1) + dx(j) )
    rslt( 3,j) = ( dx(j) + 2._r8*dx(j+1) ) / ( dx(j-1) + dx(j) )
  enddo

  !Caculate grid-based coefficients for stage 2 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    rslt( 4,j) = dx(j) / ( dx(j) + dx(j+1) )
    rslt( 5,j) = 1._r8 / sum( dx(j-1:j+2) )
    rslt( 6,j) = ( 2._r8 * dx(j+1) * dx(j) ) / ( dx(j) + dx(j+1 ) )
    rslt( 7,j) = ( dx(j-1) + dx(j  ) ) / ( 2._r8 * dx(j  ) + dx(j+1) )
    rslt( 8,j) = ( dx(j+2) + dx(j+1) ) / ( 2._r8 * dx(j+1) + dx(j  ) )
    rslt( 9,j) = dx(j  ) * ( dx(j-1) + dx(j  ) ) / ( 2._r8*dx(j  ) +    dx(j+1) )
    rslt(10,j) = dx(j+1) * ( dx(j+1) + dx(j+2) ) / (    dx(j  ) + 2._r8*dx(j+1) )
  enddo
end function compute_ppm_grids

!=======================================================================================================!



!This computes a limited parabolic interpolant using a net 5-cell stencil, but the stages of computation are broken up into 3 stages
function compute_ppm( a , dx )    result(coefs)
  !use control_mod, only: vert_remap_q_alg
  use spmd_utils            , only: masterproc
  implicit none
  real(kind=r8), intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
  real(kind=r8), intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
  real(kind=r8) ::             coefs(0:2,   nlev  )  !PPM coefficients (for parabola)
  real(kind=r8) :: ai (0:nlev  )                     !fourth-order accurate, then limited interface values
  real(kind=r8) :: dma(0:nlev+1)                     !An expression from Collela's '84 publication
  real(kind=r8) :: da                                !Ditto
  ! Hold expressions based on the grid (which are cumbersome).
  real(kind=r8) :: al, ar                            !Left and right interface values for cell-local limiting
  real(kind=r8) :: qmp,dqmax,dqmin,dq,qlc,qmin,qmax  !from Lin,2004
  real(kind=r8) :: te_top=5.6343030478959414E+05
  real(kind=r8) :: dqmono(-1:nlev+2)             !Ditto
  integer :: j
  integer :: indB, indE,vert_remap_q_alg=1,lmt
  

  ! Stage 1: Compute dma for each cell, allowing a 1-cell ghost stencil below and above the domain
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE

    da = dx(1,j) * ( dx(2,j) * ( a(j+1) - a(j) ) + dx(3,j) * ( a(j) - a(j-1) ) )
    dma(j) = minval( (/ abs(da) , 2._r8 * abs( a(j) - a(j-1) ) , 2._r8 * abs( a(j+1) - a(j) ) /) ) * sign(1._R8,da)
    if ( ( a(j+1) - a(j) ) * ( a(j) - a(j-1) ) <= 0._r8 ) dma(j) = 0._r8 

  enddo

  ! Stage 2: Compute ai for each cell interface in the physical domain (dimension nlev+1)
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    ai(j) = a(j) + dx(4,j) * ( a(j+1) - a(j) ) + dx(5,j) * ( dx(6,j) * ( dx(7,j) - dx(8,j) ) &
         * ( a(j+1) - a(j) ) - dx(9,j) * dma(j+1) + dx(10,j) * dma(j) )
  enddo

  ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain
  ! (dimension nlev) using ai on either side and ao within the cell.
  if (vert_remap_q_alg == 2) then
    indB = 3
    indE = nlev-2
  else
    indB = 1
    indE = nlev
  endif

  lmt = 1 !0 = PCM, 1 = PPM, 2 = Lin04

  if(lmt==2)then
    indB = 0
    indE = nlev+1
    !precalculating dqmono for lmt == 2 option:
    do j = indB , indE
        dq        = ( a(j+1) - a(j-1) )/4._r8
        dqmax     = max( a(j-1),a(j),a(j+1) ) - a(j)
        dqmin     = a(j) - min( a(j-1),a(j),a(j+1) )
        dqmono(j) = sign( min( abs(dqmin),dqmin,dqmax ),dq )  !mismatch (Lin,1994)
    end do

    dqmono(-1) = dqmono(0)
    dqmono(nlev+2) = dqmono(nlev+1)

    if (vert_remap_q_alg == 2) then
        indB = 3
        indE = nlev-2
    else
        indB = 1
        indE = nlev
    endif
  endif

  do j = indB , indE
    al = ai(j-1)
    ar = ai(j  )

    !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
    if(lmt==0)then
        !PCM
        coefs(0,j) = (al + ar) / 2._r8   !1.5_r8 * a(j) - ( al + ar ) / 4._r8
        coefs(1,j) = 0.0_r8 !ar - al
        coefs(2,j) = 0.0_r8 !-6._r8 * a(j) + 3._r8 * ( al + ar )
    else if (lmt==1)then
        
        !Standard PPM constraints (Colella & Woodward, 1984)    
        if ( (ar - a(j)) * (a(j) - al) <= 0._r8 ) then
                al = a(j)
                ar = a(j)
        endif
    
        if ( (ar - al) * (a(j) - (al + ar)/2._r8) >  (ar - al)**2/6._r8 ) al = 3._r8*a(j) - 2._r8 * ar
        if ( (ar - al) * (a(j) - (al + ar)/2._r8) < -(ar - al)**2/6._r8 ) ar = 3._r8*a(j) - 2._r8 * al
    
        coefs(0,j) = 1.5_r8 * a(j) - ( al + ar ) / 4._r8
        coefs(1,j) = ar - al
        coefs(2,j) = -6._r8 * a(j) + 3._r8 * ( al + ar )

        !errors are located at the start of the mapping (1 pressure) 
        !indE == top (we want bottom, indB) !problem is in indB

        !if(j==indB)then
        !        coefs(0,j) = 1.5_r8 * a(j) - ( al + ar ) / 4._r8
        !        coefs(1,j) = ar - al
        !        coefs(2,j) = -6._r8 * a(j) + 3._r8 * ( al + ar )
        !endif

!        if(j<=indB+20) then
               ! coefs(0,j) = (al + ar) / 2._r8 
               ! coefs(1,j) = 0.0_r8 
               ! coefs(2,j) = 0.0_r8 !-6._r8 * a(j) + 3._r8 * ( al + ar )
!        end if
!        if(masterproc)then
!                          if(j==1)then
!                                write (filename, '("ppm_coefs.dat")' )
!                                open(unitn, file=trim(filename),status='replace' )
!                                write(unitn,*) a(j),al,ar
!                                close(unitn)
!                          else
!                                write (filename, '("ppm_coefs.dat")' )
!                                open(unitn,file=trim(filename),status='old',position='append')
!                                write(unitn,*) a(j),al,ar
!                                close(unitn)
!                          end if
!         endif

  
    else if (lmt==2)then    

        !Lin 2004 monotonicity improvement.

        !first guess at al (q_i^-):
        
        al = ( a(j-1) + a(j) )/2._r8 + ( dqmono(j-1)-dqmono(j) )/3._r8
        ar = ai( j )        

        !Lin,2004 vertical constraint (based on Huynh,1996)

        qmp  = a(j) - 2._r8*dqmono(j)
        qlc  = a(j) + 3._r8*(dqmono(j+2) - dqmono(j))/2._r8 - dqmono(j) 
        qmin = min( a(j), qmp, qlc )
        qmax = max( a(j), qmp, qlc )
       
        al = min( max( al,qmin), qmax)
    
        qmp  = a(j) + 2._r8*dqmono(j)
        qlc  = a(j) + 3._r8*(dqmono(j) - dqmono(j-2) )/2._r8 + dqmono(j)
        qmin = min( a(j), qmp, qlc )
        qmax = max( a(j), qmp, qlc ) 
        
        ar = min( max( ar,qmin), qmax)

        coefs(2,j) = -6._r8 * a(j) + 3._r8 * ( al + ar ) 
        coefs(0,j) = 1.5_r8 * a(j) - ( al + ar ) / 4._r8 !????
        coefs(1,j) = ar - al                             !????



    end if
  enddo

  !If we're not using a mirrored boundary condition, then make the two cells bordering the top and bottom
  !material boundaries piecewise constant. Zeroing out the first and second moments, and setting the zeroth
  !moment to the cell mean is sufficient to maintain conservation.
  if (vert_remap_q_alg == 2) then
    coefs(0,1:2) = a(1:2)
    coefs(1:2,1:2) = 0._r8
    coefs(0,nlev-1:nlev) = a(nlev-1:nlev)
    coefs(1:2,nlev-1:nlev) = 0._R8
  endif
end function compute_ppm

!=======================================================================================================!


!Simple function computes the definite integral of a parabola in normalized coordinates, xi=(x-x0)/dx,
!given two bounds. Make sure this gets inlined during compilation.
function integrate_parabola( a , x1 , x2 )    result(mass)
  implicit none
  real(kind=r8), intent(in) :: a(0:2)  !Coefficients of the parabola
  real(kind=r8), intent(in) :: x1      !lower domain bound for integration
  real(kind=r8), intent(in) :: x2      !upper domain bound for integration
  real(kind=r8)             :: mass
  mass = a(0) * (x2 - x1) + a(1) * (x2 ** 2 - x1 ** 2) / 2._r8 + a(2) * (x2 ** 3 - x1 ** 3) / 3._r8
end function integrate_parabola


!=============================================================================================!

end module vertremap_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End GPU remap module    !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
