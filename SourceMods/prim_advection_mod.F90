#define _B4B 1
#define OVERLAP 1

module prim_advection_mod
!
! two formulations.  both are conservative
! u grad Q formulation:
!
!    d/dt[ Q] +  U grad Q  = 0
!
!    d/dt[ dp/dn ] = div( dp/dn U )
!
! total divergence formulation:
!    d/dt[dp/dn Q] +  div( U dp/dn Q ) = 0
!
! for convience, rewrite this as dp Q:  (since dn does not depend on time or the horizonal):
! equation is now:
!    d/dt[dp Q] +  div( U dp Q ) = 0
!
!
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use dimensions_mod,         only: nlev, np, qsize, nc
  use physconst,              only: cpair
  use derivative_mod,         only: derivative_t
  use element_mod,            only: element_t,timelevels
  use fvm_control_volume_mod, only: fvm_struct
  use hybvcoord_mod,          only: hvcoord_t
  use time_mod,               only: TimeLevel_t, TimeLevel_Qdp
  use control_mod,            only: hypervis_order, &
                                    nu_q, nu_p, limiter_option, hypervis_subcycle_q, rsplit
  use edge_mod,               only: edgevpack, edgevunpack, initedgebuffer, initedgesbuffer

  use edgetype_mod,           only: EdgeBuffer_t
  use hybrid_mod,             only: hybrid_t
  use viscosity_mod,          only: biharmonic_wk_scalar, neighbor_minmax, &
                                 neighbor_minmax_start, neighbor_minmax_finish
  use perf_mod,               only: t_startf, t_stopf, t_barrierf
  use cam_abortutils,         only: endrun
  use thread_mod,             only: horz_num_threads, tracer_num_threads

  implicit none

  private
  save

  public :: Prim_Advec_Init1, Prim_Advec_Init2
  public :: Prim_Advec_Tracers_remap
  public :: prim_advec_tracers_fvm
  public :: vertical_remap
  public :: hack_vert_rmp_init,remap_back_forth,remap_T_cons,remap_E_cons

  type (EdgeBuffer_t)      :: edgeAdv, edgeAdvp1, edgeAdvQminmax, edgeAdv1,  edgeveloc

  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  real(kind=r8), allocatable :: qmin(:,:,:), qmax(:,:,:)

!JMD I don't see why this needs to be thread private.
!JMD  type (derivative_t), public, allocatable   :: deriv(:) ! derivative struct (nthreads)
  type (derivative_t), public :: deriv


contains


  subroutine Prim_Advec_Init1(par, elem)
    use dimensions_mod, only : nlev, qsize, nelemd
    use parallel_mod,   only : parallel_t, boundaryCommMethod
    type(parallel_t)    :: par
    type (element_t)    :: elem(:)
    !
    ! Shared buffer pointers.
    ! Using "=> null()" in a subroutine is usually bad, because it makes
    ! the variable have an implicit "save", and therefore shared between
    ! threads. But in this case we want shared pointers.
    real(kind=r8), pointer :: buf_ptr(:) => null()
    real(kind=r8), pointer :: receive_ptr(:) => null()


    ! this might be called with qsize=0
    ! allocate largest one first
    ! Currently this is never freed. If it was, only this first one should
    ! be freed, as only it knows the true size of the buffer.
    call initEdgeBuffer(par,edgeAdvp1,elem,qsize*nlev + nlev,bndry_type=boundaryCommMethod,&
                         nthreads=horz_num_threads*tracer_num_threads)
    call initEdgeBuffer(par,edgeAdv,elem,qsize*nlev,bndry_type=boundaryCommMethod, &
                         nthreads=horz_num_threads*tracer_num_threads)
    call initEdgeBuffer(par,edgeAdv1,elem,nlev,bndry_type=boundaryCommMethod)
    call initEdgeBuffer(par,edgeveloc,elem,2*nlev,bndry_type=boundaryCommMethod)

    ! This is a different type of buffer pointer allocation
    ! used for determine the minimum and maximum value from
    ! neighboring  elements
    call initEdgeSBuffer(par,edgeAdvQminmax,elem,qsize*nlev*2,bndry_type=boundaryCommMethod, &
                        nthreads=horz_num_threads*tracer_num_threads)

    ! Don't actually want these saved, if this is ever called twice.
    nullify(buf_ptr)
    nullify(receive_ptr)


    ! this static array is shared by all threads, so dimension for all threads (nelemd), not nets:nete:
    allocate (qmin(nlev,qsize,nelemd))
    allocate (qmax(nlev,qsize,nelemd))

  end subroutine Prim_Advec_Init1

  subroutine Prim_Advec_Init2(fvm_corners, fvm_points)
    use dimensions_mod, only : nc
    use derivative_mod, only : derivinit

    real(kind=r8), intent(in) :: fvm_corners(nc+1)
    real(kind=r8), intent(in) :: fvm_points(nc)

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv,fvm_corners, fvm_points)
  end subroutine Prim_Advec_Init2

  !
  ! fvm driver
  !
  subroutine Prim_Advec_Tracers_fvm(elem,fvm,hvcoord,hybrid,&
        dt,tl,nets,nete)
    use fvm_consistent_se_cslam, only: run_consistent_se_cslam
    use control_mod,             only: tracer_transport_type,TRACERTRANSPORT_CONSISTENT_SE_FVM
    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (fvm_struct), intent(inout)  :: fvm(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t),        intent(in):: hybrid
    type (TimeLevel_t)                :: tl

    real(kind=r8) , intent(in) :: dt
    integer,intent(in)                :: nets,nete

    call t_barrierf('sync_prim_advec_tracers_fvm', hybrid%par%comm)
    call t_startf('prim_advec_tracers_fvm')

    if (rsplit==0) call endrun('cslam only works for rsplit>0')

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2D advection step
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (tracer_transport_type == TRACERTRANSPORT_CONSISTENT_SE_FVM) then
       call run_consistent_se_cslam(elem,fvm,hybrid,dt,tl,nets,nete,hvcoord)
    else
      call endrun('Bad tracer_transport_type in Prim_Advec_Tracers_fvm')
    end if

    call t_stopf('prim_advec_tracers_fvm')
  end subroutine Prim_Advec_Tracers_fvm



!=================================================================================================!

  subroutine Prim_Advec_Tracers_remap( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=r8) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete


    call Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
  end subroutine Prim_Advec_Tracers_remap


  subroutine euler_step_driver(np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )


  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=r8), intent(in   )         :: dt
  type (element_t)     , intent(inout)         :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier

  call euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier)

  end subroutine euler_step_driver

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! forward-in-time 2 level vertically lagrangian step
!  this code takes a lagrangian step in the horizontal
! (complete with DSS), and then applies a vertical remap
!
! This routine may use dynamics fields at timelevel np1
! In addition, other fields are required, which have to be
! explicitly saved by the dynamics:  (in elem(ie)%derived struct)
!
! Fields required from dynamics: (in
!    omega   it will be DSS'd here, for later use by CAM physics
!              we DSS omega here because it can be done for "free"
!    Consistent mass/tracer-mass advection (used if subcycling turned on)
!       dp()   dp at timelevel n0
!       vn0()  mean flux  < U dp > going from n0 to np1
!
! 3 stage
!    Euler step from t     -> t+.5
!    Euler step from t+.5  -> t+1.0
!    Euler step from t+1.0 -> t+1.5
!    u(t) = u(t)/3 + u(t+2)*2/3
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
    use derivative_mod, only : divergence_sphere
    use control_mod   , only : qsplit
    use hybrid_mod    , only : get_loop_ranges!, PrintHybrid
!    use thread_mod    , only : omp_set_num_threads, omp_get_thread_num

    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=r8) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    real (kind=r8), dimension(np,np,2     ) :: gradQ
    integer :: k,ie
    integer :: rkstage,rhs_multiplier
    integer :: n0_qdp, np1_qdp
    integer :: kbeg,kend,qbeg,qend

!    call t_barrierf('sync_prim_advec_tracers_remap_k2', hybrid%par%comm)
!    call t_startf('prim_advec_tracers_remap_rk2')
!    call extrae_user_function(1)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp) !time levels for qdp are not the same
    rkstage = 3 !   3 stage RKSSP scheme, with optimal SSP CFL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RK2 2D advection step
    ! note: stage 3 we take the oppertunity to DSS omega
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! use these for consistent advection (preserve Q=1)
    ! derived%vdp_ave        =  mean horiz. flux:   U*dp
    ! derived%omega         =  advection code will DSS this for the physics, but otherwise
    !                            it is not needed
    ! Also: save a copy of div(U dp) in derived%div(:,:,:,1), which will be DSS'd
    !       and a DSS'ed version stored in derived%div(:,:,:,2)

    call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

    do ie=nets,nete
      do k=kbeg,kend
        ! div( U dp Q),
        gradQ(:,:,1)=elem(ie)%derived%vn0(:,:,1,k)
        gradQ(:,:,2)=elem(ie)%derived%vn0(:,:,2,k)
        ! elem(ie)%derived%divdp(:,:,k) = divergence_sphere(gradQ,deriv,elem(ie))
        call divergence_sphere(gradQ,deriv,elem(ie),elem(ie)%derived%divdp(:,:,k))
        elem(ie)%derived%divdp_proj(:,:,k) = elem(ie)%derived%divdp(:,:,k)
      enddo
    enddo


    !rhs_multiplier is for obtaining dp_tracers at each stage:
    !dp_tracers(stage) = dp - rhs_multiplier*dt*divdp_proj
!    call t_startf('euler_step')

    rhs_multiplier = 0
    call euler_step_driver( np1_qdp, n0_qdp , dt/2, elem, hvcoord, hybrid, deriv, nets, nete, DSSdiv_vdp_ave, rhs_multiplier )

    rhs_multiplier = 1
    call euler_step_driver( np1_qdp, np1_qdp, dt/2, elem, hvcoord, hybrid, deriv, nets, nete, DSSno_var      , rhs_multiplier )

    rhs_multiplier = 2
    call euler_step_driver( np1_qdp, np1_qdp, dt/2, elem, hvcoord, hybrid, deriv, nets, nete, DSSomega      , rhs_multiplier )

!    call t_stopf ('euler_step')

    !to finish the 2D advection step, we need to average the t and t+2 results to get a second order estimate for t+1.
    call qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , hybrid, nets , nete )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Dissipation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( limiter_option == 8  ) then
      ! dissipation was applied in RHS.
    else
      call advance_hypervis_scalar(edgeadv,elem,hvcoord,hybrid,deriv,tl%np1,np1_qdp,nets,nete,dt)
    endif
!    call extrae_user_function(0)

!    call t_stopf('prim_advec_tracers_remap_rk2')

  end subroutine prim_advec_tracers_remap_rk2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , hybrid , nets , nete )
    use hybrid_mod, only : hybrid_t, get_loop_ranges
    implicit none
    type(element_t)     , intent(inout) :: elem(:)
    integer             , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete 
    type(hybrid_t) :: hybrid
    integer :: i,j,ie,q,k
    integer :: kbeg,kend,qbeg,qend
    real(kind=r8) :: rrkstage

    call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

#ifdef _B4B
    do ie=nets,nete
      do q=qbeg,qend
        do k=kbeg,kend
          !OMP_COLLAPSE_SIMD 
          !DIR_VECTOR_ALIGNED
          do j=1,np
          do i=1,np
             elem(ie)%state%Qdp(i,j,k,q,np1_qdp) =               &
                   ( elem(ie)%state%Qdp(i,j,k,q,n0_qdp) + &
                     (rkstage-1)*elem(ie)%state%Qdp(i,j,k,q,np1_qdp) ) / rkstage
          enddo
          enddo
        enddo
      enddo
    enddo
#else
    rrkstage=1.0_r8/real(rkstage,kind=r8)
    do ie=nets,nete
      do q=qbeg,qend
        do k=kbeg,kend
          !OMP_COLLAPSE_SIMD 
          !DIR_VECTOR_ALIGNED
          do j=1,np
          do i=1,np
           elem(ie)%state%Qdp(i,j,k,q,np1_qdp) =               &
               rrkstage *( elem(ie)%state%Qdp(i,j,k,q,n0_qdp) + &
               (rkstage-1)*elem(ie)%state%Qdp(i,j,k,q,np1_qdp) )
          enddo
          enddo
        enddo
      enddo
    enddo
#endif
  end subroutine qdp_time_avg

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS omega
  !
  ! ===================================
  use dimensions_mod , only : np, nlev
  use hybrid_mod     , only : hybrid_t!, PrintHybrid
  use hybrid_mod     , only : get_loop_ranges, threadOwnsTracer
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t, divergence_sphere, limiter_optim_iter_full
  use edge_mod       , only : edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchangev
  use hybvcoord_mod  , only : hvcoord_t

  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=r8), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier

  ! local
  real(kind=r8), dimension(np,np                       ) :: dpdiss
  real(kind=r8), dimension(np,np,nlev)                   :: dpdissk
  real(kind=r8), dimension(np,np,2                     ) :: gradQ
  real(kind=r8), dimension(np,np,2,nlev                ) :: Vstar
  real(kind=r8), dimension(np,np  ,nlev                ) :: Qtens
  real(kind=r8), dimension(np,np  ,nlev                ) :: dp
  real(kind=r8), dimension(np,np  ,nlev,qsize,nets:nete) :: Qtens_biharmonic
  real(kind=r8), dimension(np,np)                        :: div
  real(kind=r8), pointer, dimension(:,:,:)               :: DSSvar
  real(kind=r8) :: dp0(nlev)
  integer :: ie,q,i,j,k, kptr
  integer :: rhs_viss = 0
  integer :: kblk,qblk   ! The per thead size of the vertical and tracers
  integer :: kbeg, kend, qbeg, qend

  call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

  kblk = kend - kbeg + 1   ! calculate size of the block of vertical levels
  qblk = qend - qbeg + 1   ! calculate size of the block of tracers

  do k = kbeg, kend
    dp0(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
  enddo

!  call t_startf('euler_step')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss = 0
  if ( limiter_option == 8  ) then
    ! when running lim8, we also need to limit the biharmonic, so that term needs
    ! to be included in each euler step.  three possible algorithms here:
    ! 1) most expensive:
    !     compute biharmonic (which also computes qmin/qmax) during all 3 stages
    !     be sure to set rhs_viss=1
    !     cost:  3 biharmonic steps with 3 DSS
    !
    ! 2) cheapest:
    !     compute biharmonic (which also computes qmin/qmax) only on first stage
    !     be sure to set rhs_viss=3
    !     reuse qmin/qmax for all following stages (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps with 1 DSS
    !     main concern:  viscosity
    !
    ! 3)  compromise:
    !     compute biharmonic (which also computes qmin/qmax) only on last stage
    !     be sure to set rhs_viss=3
    !     compute qmin/qmax directly on first stage
    !     reuse qmin/qmax for 2nd stage stage (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps, 2 DSS
    !
    !  NOTE  when nu_p=0 (no dissipation applied in dynamics to dp equation), we should
    !        apply dissipation to Q (not Qdp) to preserve Q=1
    !        i.e.  laplace(Qdp) ~  dp0 laplace(Q)
    !        for nu_p=nu_q>0, we need to apply dissipation to Q * diffusion_dp
    !
    ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
    do ie = nets, nete
      ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
      do k = kbeg, kend
        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
        do j=1,np
        do i=1,np
          dp(i,j,k) = elem(ie)%derived%dp(i,j,k) - rhs_multiplier*dt*elem(ie)%derived%divdp_proj(i,j,k)
        enddo
        enddo
      enddo
      !JMD need to update loop based on changes in dungeon21 tag
      do q = qbeg, qend
        do k= kbeg, kend
          Qtens_biharmonic(:,:,k,q,ie) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp)/dp(:,:,k)
          if ( rhs_multiplier == 1 ) then
              qmin(k,q,ie)=min(qmin(k,q,ie),minval(Qtens_biharmonic(:,:,k,q,ie)))
              qmax(k,q,ie)=max(qmax(k,q,ie),maxval(Qtens_biharmonic(:,:,k,q,ie)))
          else
              qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
              qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
          endif
        enddo
      enddo
    enddo

    ! compute element qmin/qmax
    if ( rhs_multiplier == 0 ) then
      ! update qmin/qmax based on neighbor data for lim8
!      call t_startf('euler_neighbor_minmax1')
      call neighbor_minmax(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
!      call t_stopf('euler_neighbor_minmax1')
    endif

    ! get niew min/max values, and also compute biharmonic mixing term
    if ( rhs_multiplier == 2 ) then
      rhs_viss = 3
      ! two scalings depending on nu_p:
      ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
      ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
      if ( nu_p > 0 ) then
        do ie = nets, nete
#ifdef _B4B
          do k = kbeg, kend
            !OMP_COLLAPSE_SIMD 
            !DIR_VECTOR_ALIGNED
            do j=1,np
            do i=1,np
               dpdiss(:,:) = elem(ie)%derived%dpdiss_ave(:,:,k)
            enddo
            enddo
            do q = qbeg, qend
              ! NOTE: divide by dp0 since we multiply by dp0 below
              !OMP_COLLAPSE_SIMD
              !DIR_VECTOR_ALIGNED
              do j=1,np
              do i=1,np
                 Qtens_biharmonic(i,j,k,q,ie)=Qtens_biharmonic(i,j,k,q,ie)*dpdiss(i,j)/dp0(k)
              enddo
              enddo
            enddo
          enddo
#else
          do k = kbeg, kend
            !OMP_COLLAPSE_SIMD 
            !DIR_VECTOR_ALIGNED
            do j=1,np
            do i=1,np
               dpdissk(i,j,k) = elem(ie)%derived%dpdiss_ave(i,j,k)/dp0(k)
            enddo
            enddo
          enddo
          do q = qbeg,qend
            do k = kbeg, kend
              ! NOTE: divide by dp0 since we multiply by dp0 below
              !OMP_COLLAPSE_SIMD 
              !DIR_VECTOR_ALIGNED
              do j=1,np
              do i=1,np
                 Qtens_biharmonic(i,j,k,q,ie)=Qtens_biharmonic(i,j,k,q,ie)*dpdissk(i,j,k)
              enddo
              enddo
            enddo
          enddo
#endif
        enddo
      endif

!   Previous version of biharmonic_wk_scalar_minmax included a min/max
!   calculation into the boundary exchange.  This was causing cache issues.
!   Split the single operation into two separate calls
!      call neighbor_minmax()
!      call biharmonic_wk_scalar()
!
!      call biharmonic_wk_scalar_minmax( elem , qtens_biharmonic , deriv , edgeAdvQ3 , hybrid , &
!           nets , nete , qmin(:,:,nets:nete) , qmax(:,:,nets:nete) )
#ifdef OVERLAP
      call neighbor_minmax_start(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
      call biharmonic_wk_scalar(elem,qtens_biharmonic,deriv,edgeAdv,hybrid,nets,nete)
      do ie = nets, nete
        do q = qbeg, qend
          do k = kbeg, kend
            !OMP_COLLAPSE_SIMD 
            !DIR_VECTOR_ALIGNED
            do j=1,np
            do i=1,np
               ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
               qtens_biharmonic(i,j,k,q,ie) = &
                     -rhs_viss*dt*nu_q*dp0(k)*Qtens_biharmonic(i,j,k,q,ie) / elem(ie)%spheremp(i,j)
            enddo
            enddo
          enddo
        enddo
      enddo
      call neighbor_minmax_finish(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
#else
      call t_startf('euler_neighbor_minmax2')
      call neighbor_minmax(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
      call t_stopf('euler_neighbor_minmax2')
      call biharmonic_wk_scalar(elem,qtens_biharmonic,deriv,edgeAdv,hybrid,nets,nete)

      do ie = nets, nete
        do q = qbeg, qend
          do k = kbeg, kend
            !OMP_COLLAPSE_SIMD 
            !DIR_VECTOR_ALIGNED
            do j=1,np
            do i=1,np
               ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
               qtens_biharmonic(i,j,k,q,ie) = &
                     -rhs_viss*dt*nu_q*dp0(k)*Qtens_biharmonic(i,j,k,q,ie) / elem(ie)%spheremp(i,j)
            enddo
            enddo
          enddo
        enddo
      enddo
#endif


    endif
  endif  ! compute biharmonic mixing term and qmin/qmax


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie = nets, nete


    ! Compute velocity used to advance Qdp
    do k = kbeg, kend
      ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
      ! but that's ok because rhs_multiplier=0 on the first stage:
      !OMP_COLLAPSE_SIMD 
      !DIR_VECTOR_ALIGNED
      do j=1,np
      do i=1,np
         dp(i,j,k) = elem(ie)%derived%dp(i,j,k) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(i,j,k)
         Vstar(i,j,1,k) = elem(ie)%derived%vn0(i,j,1,k) / dp(i,j,k)
         Vstar(i,j,2,k) = elem(ie)%derived%vn0(i,j,2,k) / dp(i,j,k)
      enddo
      enddo
    enddo
    if ( limiter_option == 8) then
        ! Note that the term dpdissk is independent of Q
        do k = kbeg, kend
          ! UN-DSS'ed dp at timelevel n0+1:
          !OMP_COLLAPSE_SIMD 
          !DIR_VECTOR_ALIGNED
          do j=1,np
          do i=1,np
             dpdissk(i,j,k) = dp(i,j,k) - dt * elem(ie)%derived%divdp(i,j,k)
          enddo
          enddo
          if ( nu_p > 0 .and. rhs_viss /= 0 ) then
            ! add contribution from UN-DSS'ed PS dissipation
!            dpdiss(:,:) = ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) *
!            elem(ie)%derived%psdiss_biharmonic(:,:)
            !OMP_COLLAPSE_SIMD 
            !DIR_VECTOR_ALIGNED
            do j=1,np
            do i=1,np
               dpdiss(i,j) = elem(ie)%derived%dpdiss_biharmonic(i,j,k)
               dpdissk(i,j,k) = dpdissk(i,j,k) - rhs_viss * dt * nu_q * dpdiss(i,j) / elem(ie)%spheremp(i,j)
            enddo
            enddo
          endif
          ! IMPOSE ZERO THRESHOLD.  do this here so it can be turned off for
          ! testing
          do q=qbeg, qend
             qmin(k,q,ie)=max(qmin(k,q,ie),0.0_r8)
          enddo
        enddo
    endif  ! limiter == 8


    ! advance Qdp
    do q = qbeg, qend
      do k = kbeg, kend
        ! div( U dp Q),
        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
        do j=1,np
        do i=1,np
           gradQ(i,j,1) = Vstar(i,j,1,k) * elem(ie)%state%Qdp(i,j,k,q,n0_qdp)
           gradQ(i,j,2) = Vstar(i,j,2,k) * elem(ie)%state%Qdp(i,j,k,q,n0_qdp)
        enddo
        enddo
        ! Qtens(:,:,k) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp) - &
        !               dt * divergence_sphere( gradQ , deriv , elem(ie) )
        call divergence_sphere( gradQ , deriv , elem(ie),div )

        !OMP_COLLAPSE_SIMD
        !DIR_VECTOR_ALIGNED
        do j=1,np
          do i=1,np
            Qtens(i,j,k) = elem(ie)%state%Qdp(i,j,k,q,n0_qdp) - dt * div(i,j)
          enddo
        enddo

        ! optionally add in hyperviscosity computed above:
        if ( rhs_viss /= 0 ) then 
          !OMP_COLLAPSE_SIMD 
          !DIR_VECTOR_ALIGNED
          do j=1,np
          do i=1,np
              Qtens(i,j,k) = Qtens(i,j,k) + Qtens_biharmonic(i,j,k,q,ie)
          enddo
          enddo
        endif
      enddo

      if ( limiter_option == 8) then
        ! apply limiter to Q = Qtens / dp_star
        call limiter_optim_iter_full( Qtens(:,:,:) , elem(ie)%spheremp(:,:) , qmin(:,q,ie) , &
                                                          qmax(:,q,ie) , dpdissk, kbeg, kend )
      endif


      ! apply mass matrix, overwrite np1 with solution:
      ! dont do this earlier, since we allow np1_qdp == n0_qdp
      ! and we dont want to overwrite n0_qdp until we are done using it
      do k = kbeg, kend
        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
        do j=1,np
        do i=1,np
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%spheremp(i,j) * Qtens(i,j,k)
        enddo
        enddo
      enddo

      if ( limiter_option == 4 ) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        ! sign-preserving limiter, applied after mass matrix
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!JMD !$OMP BARRIER
!JMD !$OMP MASTER
        call limiter2d_zero(elem(ie)%state%Qdp(:,:,:,q,np1_qdp))
!JMD !$OMP END MASTER
!JMD !$OMP BARRIER
      endif

      kptr = nlev*(q-1) + kbeg - 1
      call edgeVpack(edgeAdvp1 , elem(ie)%state%Qdp(:,:,kbeg:kend,q,np1_qdp) , kblk , kptr , ie )
    enddo
    ! only perform this operation on thread which owns the first tracer
    if (DSSopt>0) then
      if (threadOwnsTracer(hybrid,1)) then
        ! all zero so we only have to DSS 1:nlev
        if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega(:,:,:)
        if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
        ! also DSS extra field
        do k = kbeg, kend
          !OMP_COLLAPSE_SIMD 
          !DIR_VECTOR_ALIGNED
          do j=1,np
          do i=1,np
             DSSvar(i,j,k) = elem(ie)%spheremp(i,j) * DSSvar(i,j,k)
          enddo
          enddo
        enddo
        kptr = nlev*qsize + kbeg - 1
        call edgeVpack( edgeAdvp1 , DSSvar(:,:,kbeg:kend), kblk, kptr, ie)
      endif
    end if
  enddo

  call bndry_exchangeV( hybrid , edgeAdvp1)

  do ie = nets, nete
    ! only perform this operation on thread which owns the first tracer
    if (DSSopt>0) then
      if(threadOwnsTracer(hybrid,1)) then
        if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega(:,:,:)
        if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
        kptr = qsize*nlev + kbeg -1
        call edgeVunpack( edgeAdvp1 , DSSvar(:,:,kbeg:kend) , kblk , kptr , ie )
        do k = kbeg, kend
           !OMP_COLLAPSE_SIMD 
           !DIR_VECTOR_ALIGNED
           do j=1,np
           do i=1,np
               DSSvar(i,j,k) = DSSvar(i,j,k) * elem(ie)%rspheremp(i,j)
           enddo
           enddo
        enddo
      endif
    end if
    do q = qbeg, qend
      kptr = nlev*(q-1) + kbeg - 1
      call edgeVunpack( edgeAdvp1 , elem(ie)%state%Qdp(:,:,kbeg:kend,q,np1_qdp) , kblk , kptr , ie )
        do k = kbeg, kend
          !OMP_COLLAPSE_SIMD 
          !DIR_VECTOR_ALIGNED
          do j=1,np
          do i=1,np
             elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%rspheremp(i,j) * elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
          enddo
          enddo
        enddo
    enddo
  enddo
!  call t_stopf('euler_step')

  end subroutine euler_step



  subroutine limiter2d_zero(Q)
  ! mass conserving zero limiter (2D only).  to be called just before DSS
  !
  ! this routine is called inside a DSS loop, and so Q had already
  ! been multiplied by the mass matrix.  Thus dont include the mass
  ! matrix when computing the mass = integral of Q over the element
  !
  ! ps is only used when advecting Q instead of Qdp
  ! so ps should be at one timelevel behind Q
  implicit none
  real (kind=r8), intent(inout) :: Q(np,np,nlev)

  ! local
!  real (kind=r8) :: dp(np,np)
  real (kind=r8) :: mass,mass_new,ml
  integer i,j,k

  do k = nlev , 1 , -1
    mass = 0
    do j = 1 , np
      do i = 1 , np
        !ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
        ml = Q(i,j,k)
        mass = mass + ml
      enddo
    enddo

    ! negative mass.  so reduce all postive values to zero
    ! then increase negative values as much as possible
    if ( mass < 0 ) Q(:,:,k) = -Q(:,:,k)
    mass_new = 0
    do j = 1 , np
      do i = 1 , np
        if ( Q(i,j,k) < 0 ) then
          Q(i,j,k) = 0
        else
          ml = Q(i,j,k)
          mass_new = mass_new + ml
        endif
      enddo
    enddo

    ! now scale the all positive values to restore mass
    if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass) / mass_new
    if ( mass     < 0 ) Q(:,:,k) = -Q(:,:,k)
  enddo
  end subroutine limiter2d_zero

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine advance_hypervis_scalar( edgeAdv , elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
  !  hyperviscsoity operator for foward-in-time scheme
  !  take one timestep of:
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  use dimensions_mod , only : np, nlev
  use hybrid_mod     , only : hybrid_t!, PrintHybrid
  use hybrid_mod     , only : get_loop_ranges
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t
  use edge_mod       , only : edgevpack, edgevunpack
  use edgetype_mod   , only : EdgeBuffer_t
  use bndry_mod      , only : bndry_exchangev

  implicit none
  type (EdgeBuffer_t)  , intent(inout)         :: edgeAdv
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nt
  integer              , intent(in   )         :: nt_qdp
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  real (kind=r8), intent(in   )         :: dt2

  ! local
  real (kind=r8), dimension(np,np,nlev,qsize,nets:nete) :: Qtens
  real (kind=r8), dimension(np,np,nlev                ) :: dp
!  real (kind=r8), dimension(      nlev,qsize,nets:nete) :: min_neigh
!  real (kind=r8), dimension(      nlev,qsize,nets:nete) :: max_neigh
  integer :: k,kptr,ie,ic,q,i,j
  integer :: kbeg,kend,qbeg,qend

! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
!       data is incorrect (offset by a few numbers actually)
!       removed for now.
!  real (kind=r8), dimension(:,:), pointer :: spheremp,rspheremp
  real (kind=r8) :: dt,dp0
  integer :: kblk,qblk   ! The per thead size of the vertical and tracers

  call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

  if ( nu_q           == 0 ) return
  if ( hypervis_order /= 2 ) return

  kblk = kend - kbeg + 1   ! calculate size of the block of vertical levels
  qblk = qend - qbeg + 1   ! calculate size of the block of tracers

  call t_startf('advance_hypervis_scalar')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dt = dt2 / hypervis_subcycle_q

  do ic = 1 , hypervis_subcycle_q
    do ie = nets, nete
      ! Qtens = Q/dp   (apply hyperviscsoity to dp0 * Q, not Qdp)
      do k = kbeg, kend
         ! various options:
         !   1)  biharmonic( Qdp )
         !   2)  dp0 * biharmonic( Qdp/dp )
         !   3)  dpave * biharmonic(Q/dp)
         ! For trace mass / mass consistenciy, we use #2 when nu_p=0
         ! and #e when nu_p>0, where dpave is the mean mass flux from the nu_p
         ! contribution from dynamics.
         dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * hvcoord%ps0
         dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
         if (nu_p>0) then
            do q = qbeg, qend
               Qtens(:,:,k,q,ie) = elem(ie)%derived%dpdiss_ave(:,:,k)*&
                    elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
            enddo
         else
            do q = qbeg, qend
               Qtens(:,:,k,q,ie) = dp0*elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
            enddo
         endif
      enddo
   enddo

    ! compute biharmonic operator. Qtens = input and output
    call biharmonic_wk_scalar( elem , Qtens , deriv , edgeAdv , hybrid ,  nets , nete )

    do ie = nets, nete
      !spheremp     => elem(ie)%spheremp
      do q = qbeg, qend
        do k = kbeg, kend
          dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * hvcoord%ps0
         do j = 1 , np
            do i = 1 , np

              ! advection Qdp.  For mass advection consistency:
              ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )
              elem(ie)%state%Qdp(i,j,k,q,nt_qdp) = elem(ie)%state%Qdp(i,j,k,q,nt_qdp) * elem(ie)%spheremp(i,j) &
                                                                                 - dt * nu_q * Qtens(i,j,k,q,ie)
        enddo
        enddo
        enddo

        if (limiter_option .ne. 0 ) then
!JMD Only need if threading over the vertical
!JMD!$OMP BARRIER
!JMD!$OMP MASTER
           ! smooth some of the negativities introduced by diffusion:
           call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,nt_qdp) )
!JMD!$OMP END MASTER
!JMD!$OMP BARRIER
        endif

      enddo
      do q = qbeg, qend
         kptr = nlev*(q-1) + kbeg - 1
         call edgeVpack( edgeAdv , elem(ie)%state%Qdp(:,:,kbeg:kend,q,nt_qdp) , kblk, kptr, ie )
      enddo
    enddo

    call bndry_exchangeV( hybrid , edgeAdv)

    do ie = nets, nete
      do q = qbeg, qend
         kptr = nlev*(q-1) + kbeg - 1
         call edgeVunpack( edgeAdv , elem(ie)%state%Qdp(:,:,kbeg:kend,q,nt_qdp) , kblk, kptr, ie )
      enddo
      !rspheremp     => elem(ie)%rspheremp
      do q = qbeg, qend
        ! apply inverse mass matrix
        do k = kbeg, kend
          elem(ie)%state%Qdp(:,:,k,q,nt_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,nt_qdp)
        enddo
      enddo
    enddo

  enddo
  call t_stopf('advance_hypervis_scalar')
  end subroutine advance_hypervis_scalar


  subroutine vertical_remap(hybrid,elem,fvm,hvcoord,dt,np1,np1_qdp,n_fvm,nets,nete)
    ! This routine is called at the end of the vertically Lagrangian
    ! dynamics step to compute the vertical flux needed to get back
    ! to reference eta levels
    !
    ! input:
    !     derived%dp()  delta p on levels at beginning of timestep
    !     state%dp3d(np1)  delta p on levels at end of timestep
    ! output:
    !     state%psdry(np1)          surface pressure at time np1
    !i
    
    use hybvcoord_mod, only          : hvcoord_t
    use vertremap_mod,          only : remap1, remap1_nofilter,remap_Q_ppm
    use hybrid_mod            , only : hybrid_t!, set_region_num_threads
    use fvm_control_volume_mod, only : fvm_struct
    use control_mod,            only : se_prescribed_wind_2d
    use dimensions_mod        , only : ntrac
    use dimensions_mod        , only : qsize_condensate_loading, qsize_condensate_loading_idx
    use hack_vert_rmp        ! , only : get_levels,overwrite_state,write_data,write_data_TE,diagnostic,lhack_vert_rmp,diagnostic_eul,write_data_eul
    use spmd_utils            , only: masterproc
    !use physconst,              only: r_universal
    use MOM_remapping

    type (hybrid_t),  intent(in)    :: hybrid  ! distributed parallel structure (shared)
    type(fvm_struct), intent(inout) :: fvm(:)
    type (element_t), intent(inout) :: elem(:)
    integer,          intent(in)    :: n_fvm
    !
    real (kind=r8) :: cdp(1:nc,1:nc,nlev,ntrac)
    real (kind=r8) :: dpc(nc,nc,nlev),dpc_star(nc,nc,nlev)
    
    type (hvcoord_t) :: hvcoord
    real (kind=r8)   :: dt
    integer          :: ie,i,j,k,np1,nets,nete,np1_qdp,q, m_cnst,temp_int,unitn
    character(len=256):: filename
    real (kind=r8), dimension(np,np,nlev)  :: dp_moist,dp_star_moist, dp_inv,dp_dry,dp_star_dry,dp_s_inv
    real (kind=r8), dimension(np,np,nlev)  :: q_test,q_test_s
    real (kind=r8), dimension(np,np,nlev,2):: ttmp
    real (kind=r8), dimension(np,np,nlev+1):: phi
    real(KIND=r8), dimension(nlev+1) :: pint1,pint2,dx1
    real(KIND=r8), dimension(nlev) :: t0,lnp1,lnp2
    logical          :: remap_te,filtered,phi_inc,ppm,pqm  
    real(KIND=r8) :: r_universal
    real(KIND=r8) :: dpphi,kappa
    real(KIND=r8),dimension(np,np,nlev) :: num_tmp,den_tmp
    type(remapping_CS)                  :: CS
    real, dimension(nlev) :: dp_mom,ttmp_mom,dp_s_mom,ttmp_s_mom,qdp_tmp,qdp_s_tmp


    real (kind=r8), dimension(np,np,nlev,2) :: E_0
    real (kind=r8), dimension(np,np,nlev,2) :: E_1

    real (kind=r8), dimension(np,np,2,nlev,3) :: v_0 
    real (kind=r8), dimension(np,np,2,nlev,3) :: v_1

    real (kind=r8), dimension(np,np,nlev,3) :: t_0  
    real (kind=r8), dimension(np,np,nlev,3) :: t_1


    ! reference levels:
    !   dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps(i,j)
    !   hybi(1)=0          pure pressure at top of atmosphere
    !   hyai(1)=ptop
    !   hyai(nlev+1) = 0   pure sigma at bottom
    !   hybi(nlev+1) = 1
    !
    ! sum over k=1,nlev
    !  sum(dp(k)) = (hyai(nlev+1)-hyai(1))*ps0 + (hybi(nlev+1)-hybi(1))*ps
    !             = -ps0 + ps
    !  ps =  ps0+sum(dp(k))
    !
    ! reference levels:
    !    dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps
    !
    do ie=nets,nete
      !
      !  REMAP u,v,T from levels in dp3d() to REF levels
      !
      ! update final ps
      elem(ie)%state%psdry(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
           sum(elem(ie)%state%dp3d(:,:,:,np1),3)
      !

      if(lhack_vert_rmp) then                           !Remapping fixed profile
           
        ttmp = 0._r8
        phi = 0._r8
        dpphi = 0._r8
        q_test = 1.0_r8
        q_test(:,:,15) = 2.0_r8
        r_universal = 1.38065e-23*6.02214e26/28.966
        kappa = r_universal/cpair

        call hack_vert_rmp_init(dp_star_moist,dp_moist,elem(ie)%state%t,elem(ie)%state%v,pint1,pint2,lnp1,lnp2,1,1)

        if (minval(dp_star_moist)<0) call endrun('negative moist layer thickness.  timestep or remap time too large')
        if (minval(dp_moist)<0) call endrun('negative layer thickness.  timestep or remap time too large')

        dp_inv=1.0_R8/dp_moist !for efficiency
        dp_s_inv=1.0_R8/dp_star_moist

        call write_data(pint1,elem(ie)%state%t(1,1,:,np1),elem(ie)%state%v(1,1,1,:,np1),0) !Print state 0

        call diagnostic(dp_star_moist*elem(ie)%state%v(:,:,1,:,np1),&
                dp_star_moist*elem(ie)%state%t(:,:,:,np1)*cpair,&
                dp_star_moist*((elem(ie)%state%v(:,:,1,:,np1)**2+elem(ie)%state%v(:,:,2,:,np1)**2)/2._r8),dpphi,0)  !diagnostic stage 0 
                
                !END hacked-profile initialization

      else      !Preserving old code
        do k=1,nlev
                dp_dry(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                        ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%psdry(:,:,np1)
                dp_star_dry(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
        enddo
        if (minval(dp_star_dry)<0) call endrun('negative dry layer thickness.  timestep or remap time too large B')
      !     
        dp_star_moist(:,:,:) = dp_star_dry(:,:,:)
        do q=1,qsize_condensate_loading
                m_cnst = qsize_condensate_loading_idx(q)
                do k=1,nlev
                        dp_star_moist(:,:,k)= dp_star_moist(:,:,k)+elem(ie)%state%Qdp(:,:,k,m_cnst,np1_qdp)
                end do
        end do
        if (minval(dp_star_moist)<0) call endrun('negative moist layer thickness.  timestep or remap time too large')
 
        call remap1(elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp),np,1,qsize,qsize,dp_star_dry,dp_dry,hybrid=hybrid)
  
        !
        ! compute moist reference pressure level thickness
       
        dp_moist(:,:,:) = dp_dry(:,:,:)
        do q=1,qsize_condensate_loading
                m_cnst = qsize_condensate_loading_idx(q)
                do k=1,nlev
                        dp_moist(:,:,k) = dp_moist(:,:,k)+elem(ie)%state%Qdp(:,:,k,m_cnst,np1_qdp)
                end do
        end do
        if (minval(dp_star_moist)<0) call endrun('negative layer thickness.  timestep or remap time too large')
        !Possible bug: should it be dp_moist here?
        dp_inv=1.0_R8/dp_moist !for efficiency
      end if                

      !!TRIVIAL CASE -- SANITY CHECK:
!      dp_moist = dp_star_moist
!      lnp2 = lnp1
!      dp_inv = dp_s_inv


      remap_te = .true.   !using Total Energy remapping instead of Temperature    
      !filtered = .false.   !using remap1() instead of remap1_nofilter()
      !ppm = .false.        !CAM PPM
      pqm = .true.       !MOM REMAPS
      !phi_inc  = .true.   !using Geopotential term inside TE

      if(pqm)then
           call initialize_remapping(CS,'PPM_H4',.true.,.true.,.true.,.true.)
           if(masterproc)then
             print *, "INIT MOM REMAP PASSED..."
           endif
      endif

      if(lhack_vert_rmp) then

!      do k=1,nlev
!           dp_dry(:,:,k) = dp_moist(:,:,k)
!           dp_star_dry(:,:,k) = dp_star_moist(:,:,k)
!      enddo

      !qdp_s_tmp = q_test(1,1,:)               



          if(remap_te)then
call remap_back_forth(99,1,3,dp_star_moist,ttmp,elem(ie)%state%t,elem(ie)%state%v,pint1,dp_moist,E_1,t_1,v_1,pint2,elem,np1)                              
          else!if(.not.remap_te)then
call remap_back_forth(99,2,3,dp_star_moist,ttmp,elem(ie)%state%t,elem(ie)%state%v,pint1,dp_moist,E_1,t_1,v_1,pint2,elem,np1)
          end if
      

          exit

   !if not hack-remapping:
!*********************************************************************************************************
      else
          !Need to compute the interfaces pint1,2:

          do k=1,nlev+1
                pint1(k) = 0.00000000001_r8
                pint2(k) = 0.00000000001_r8
                do i = 1,k-1
                        pint1(k) = pint1(k) + dp_star_moist(1,1,i)
                        pint2(k) = pint2(k) + dp_moist(1,1,i)
                enddo

                if(masterproc)then
                      print *, 'PINT LEVELS :',pint1(k),pint2(k),k
                      write (filename, '("stat_", I0.3, ".dat")' )  k
                      open(unitn, file=trim(filename), status='replace' )
                          write(unitn,*) pint1(k),pint2(k),k
                      close(unitn)
                endif
          enddo          
   
          if(remap_te)then
                ! remap u,v and cpair*T + .5 u^2
!                ttmp(:,:,:,1)=(elem(ie)%state%v(:,:,1,:,np1)**2 + &
!                        elem(ie)%state%v(:,:,2,:,np1)**2)/2 + &
!                        elem(ie)%state%t(:,:,:,np1)*cpair
!          ttmp = 1._r8
!          pint1 = 1._r8
!          pint2 = 2._r8
call remap_E_cons(1,dp_star_moist,ttmp,elem(ie)%state%t,elem(ie)%state%v,pint1,dp_moist,E_1,t_1,v_1,pint2,elem,np1)
          ttmp = E_1
          elem(ie)%state%t = t_1
          elem(ie)%state%v = v_1
          else

call remap_T_cons(2,dp_star_moist,elem(ie)%state%t,elem(ie)%state%v,pint1,dp_moist,t_1,v_1,pint2,elem,np1)
          elem(ie)%state%t = t_1
          elem(ie)%state%v = v_1

!                ttmp(:,:,:,1)=elem(ie)%state%t(:,:,:,np1)
          end if

!          ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_star_moist
!          call remap1(ttmp,np,1,1,1,dp_star_moist,dp_moist)
!          elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)*dp_inv  
      
!          ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star_moist  !u
!          ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star_moist  !v
!          call remap1(ttmp,np,1,2,2,dp_star_moist,dp_moist)
          !        call remap1_nofilter(ttmp,np,2,dp_star,dp)

!          if ( .not. se_prescribed_wind_2d ) &
!                 elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)*dp_inv    !u_rmp
!          if ( .not. se_prescribed_wind_2d ) &
!                 elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)*dp_inv    !v_rmp
     
!          if(remap_te)then
                ! back out T from TE
!                elem(ie)%state%t(:,:,:,np1) = &
!                ( elem(ie)%state%t(:,:,:,np1) - ((elem(ie)%state%v(:,:,1,:,np1)**2 + &
!                elem(ie)%state%v(:,:,2,:,np1)**2)/2))/cpair           !T_rmp
!          end if

      end if 


      ! remap the gll tracers from lagrangian levels (dp_star)  to REF levels dp
      if (qsize>0) then
        
        if ( se_prescribed_wind_2d ) then
          ! Peter Lauritzen et al, "The terminator 'toy'-chemistry test: A simple tool to assess errors in transport schemes",
          !   submitted to Geosci Model Dev, Oct 2014
          ! -- code to let dp evolve without vertical transport of tracers (consistent mass tracer coupling)
          do q=1,qsize
            do k=1,nlev
              do j=1,np
                do i=1,np
                  !elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp) * dp(i,j,k)/dp_star(i,j,k)
                  ttmp(i,j,k,1)= elem(ie)%state%Qdp(i,j,k,q,np1_qdp) / dp_star_moist(i,j,k) ! This is the actual q
                  elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = ttmp(i,j,k,1) * dp_moist(i,j,k)
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      
      
      if (ntrac>0) then
        do i=1,nc
          do j=1,nc
            !
            ! compute source (cdp) and target (dpc) pressure grids for vertical remapping
            !
            do k=1,nlev
              dpc(i,j,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                   (hvcoord%hybi(k+1) - hvcoord%hybi(k))*fvm(ie)%psc(i,j)
              cdp(i,j,k,1:ntrac)=fvm(ie)%c(i,j,k,1:ntrac,n_fvm)*fvm(ie)%dp_fvm(i,j,k,n_fvm)
            end do
          end do
        end do
        dpc_star=fvm(ie)%dp_fvm(1:nc,1:nc,:,n_fvm)
        call remap1(cdp,nc,1,ntrac,ntrac,dpc_star,dpc)
        do k=1,nlev
          do j=1,nc
            do i=1,nc
              fvm(ie)%dp_fvm(i,j,k,n_fvm)=dpc(i,j,k)
              fvm(ie)%c(i,j,k,1:ntrac,n_fvm)=cdp(i,j,k,1:ntrac)/dpc(i,j,k)
            end do
          end do
        end do
      end if
      !         call remap_velocityC(np1,dt,elem,fvm,hvcoord,ie)
    enddo

  end subroutine vertical_remap

  subroutine hack_vert_rmp_init(dp_star_moist,dp_moist,t,v,pint1,pint2,lnp1,lnp2,num,vel_0)

  use shr_kind_mod,           only: r8=>shr_kind_r8
  use hack_vert_rmp  

    use hybvcoord_mod, only          : hvcoord_t
    use vertremap_mod,          only : remap1, remap1_nofilter,remap_Q_ppm
    use hybrid_mod            , only : hybrid_t!, set_region_num_threads
    use fvm_control_volume_mod, only : fvm_struct
    use control_mod,            only : se_prescribed_wind_2d
    use dimensions_mod        , only : ntrac
    use dimensions_mod        , only : qsize_condensate_loading,qsize_condensate_loading_idx
    use hack_vert_rmp
    use spmd_utils            , only: masterproc
 

  real (kind=r8), dimension(np,np,nlev), intent(inout)  :: dp_moist,dp_star_moist
  real (kind=r8), dimension(nlev), intent(out) :: lnp1,lnp2 
  real (kind=r8), dimension(nlev) :: t0
  real (kind=r8), dimension(np,np,nlev,3),intent(inout) :: t
  real (kind=r8), dimension(np,np,2,nlev,3),intent(inout) :: v
  integer          :: i,j,k
  integer, intent(in) :: num,vel_0
  real(KIND=r8) :: dpphi,kappa
  real(KIND=r8), dimension(nlev+1), intent(inout) :: pint1,pint2

        call get_levels(num,pint1,pint2,t0) !option 1 is CAM5

        dp_star_moist = 0._r8 !initializing as zeros
        dp_moist = 0._r8

        do k=1,size(t0)
                 
                 dp_star_moist(:,:,k) = pint1(k+1)-pint1(k)
                 lnp1(k) = log(pint1(k+1))-log(pint1(k))
        
                 dp_moist(:,:,k) = pint2(k+1)-pint2(k)
                 lnp2(k) = log(pint2(k+1))-log(pint2(k))

                 t(:,:,k,:) = t0(k)
        end do

        if(vel_0==1)then
                !setting velocities to zero:
                 v(:,:,1,:,:) = 0._r8
                 v(:,:,2,:,:) = 0._r8
        endif

end subroutine hack_vert_rmp_init

  subroutine remap_E_cons(rmp_kind,dp_0,E_0,t_0,v_0,pint_0,dp_1,E_1,t_1,v_1,pint_1,elem,np1)
  !This routine will remap an energy field from system 0 to system 1
  !conserving it's total budget, with different remapping algorithms,
  !following the total energy formula from Lin'04

  !I will add all the modules i was using in vertical_remap hoping for the best:

    use hybvcoord_mod, only          : hvcoord_t
    use vertremap_mod,          only : remap1, remap1_nofilter,remap_Q_ppm
    use hybrid_mod            , only : hybrid_t!, set_region_num_threads
    use fvm_control_volume_mod, only : fvm_struct
    use control_mod,            only : se_prescribed_wind_2d
    use dimensions_mod        , only : ntrac
    use dimensions_mod        , only : qsize_condensate_loading,qsize_condensate_loading_idx
    use hack_vert_rmp        
    use spmd_utils            , only: masterproc
    !use physconst,              only: r_universal
    use MOM_remapping

    type (element_t), intent(inout) :: elem(:)

   !
    type (hvcoord_t) :: hvcoord
    real (kind=r8)   :: dt
    integer          :: ie,i,j,k,nets,nete,np1_qdp,q, m_cnst,temp_int,unitn
    character(len=256):: filename
    real (kind=r8), dimension(np,np,nlev)  :: dp_1_inv,dp_0_inv
    real (kind=r8), dimension(np,np,nlev)  :: q_test,q_test_s
    real (kind=r8), dimension(np,np,nlev,2):: ttmp
    real (kind=r8), dimension(np,np,nlev+1):: phi
    real(KIND=r8), dimension(nlev+1) :: pint1,pint2,dx1
    real(KIND=r8), dimension(nlev) :: t0,lnp1,lnp2
    logical          :: remap_te,filtered,phi_inc,ppm,pqm
    real(KIND=r8) :: r_universal
    real(KIND=r8) :: dpphi,kappa
    real(KIND=r8),dimension(np,np,nlev) :: num_tmp,den_tmp
    type(remapping_CS)                  :: CS
    real, dimension(nlev) :: dp_mom,ttmp_mom,dp_s_mom,ttmp_s_mom,qdp_tmp,qdp_s_tmp

    integer, intent(in) :: rmp_kind,np1

    real (kind=r8), dimension(np,np,nlev,2),intent(in) :: E_0
    real (kind=r8), dimension(np,np,nlev,2),intent(inout) :: E_1

    real (kind=r8), dimension(np,np,2,nlev,3),intent(in) :: v_0  !Vel. last arg is 'timelevels'
    real (kind=r8), dimension(np,np,2,nlev,3),intent(out) :: v_1

    real (kind=r8), dimension(np,np,nlev,3),intent(in) :: t_0  !Temp. last arg same
    real (kind=r8), dimension(np,np,nlev,3),intent(out) :: t_1

    real (kind=r8), dimension(np,np,nlev),intent(in) :: dp_0,dp_1
    real (kind=r8), dimension(nlev+1),intent(in) :: pint_0,pint_1



    !============
!    rmp_kind(s):
    ! 1.- filtered (remap1_filter())
    ! 2.- ppm (remap_Q_ppm)
    ! 3.- mom (remapping_core_h() using method from vertical_remap())
    ! 4.- no filter (remap1_nofilter())
    !============
 
    unitn = 8 
    r_universal = 1.38065e-23*6.02214e26/28.966
    kappa = r_universal/cpair

    phi_inc = .true.
    phi = 0._r8
    
    dp_0_inv = 1._r8/dp_0
    dp_1_inv = 1._r8/dp_1

    !Calculating energy at system 0
    ttmp(:,:,:,1) = (v_0(:,:,1,:,np1)**2 + v_0(:,:,2,:,np1)**2)/2._r8 + t_0(:,:,:,np1)*cpair 

    if(phi_inc)then
    !Calculating Phi term:
    do k=nlev,1,-1

        if(masterproc)then
              write(*,*) "DEBUG LEVELS :", phi(1,1,k),ttmp(1,1,k,1),pint_0(k+1),pint_0(k),k
              !write (filename, '("stat_0_", I0.3, ".dat")' )  k
              !open(unitn, file=trim(filename), status='replace' )
                ! write(unitn,*) phi(1,1,k),ttmp(1,1,k,1),pint_0(k+1),pint_0(k)!
              
              !close(unitn)
        endif

        phi(:,:,k) = phi(:,:,k+1) + r_universal*t_0(:,:,k,np1)*( log(pint_0(k+1)) - log(pint_0(k)) )

        ttmp(:,:,k,1)= ttmp(:,:,k,1) + (pint_0(k+1)*phi(:,:,k+1) - &
                       pint_0(k)*phi(:,:,k))*dp_0_inv(:,:,k) !Energy  

    end do
    end if

    ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_0 !E*dp_star_moist 

    if(rmp_kind==1)then
         call remap1(ttmp,np,1,1,1,dp_0,dp_1) !E_rmp*dp_moist
    elseif(rmp_kind==2)then
         call remap_Q_ppm(ttmp,np,1,1,1,dp_0,dp_1,1)
!         call remap_Q_ppm(elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp),np,1,qsize,qsize,dp_star_dry,dp_dry,2)
    elseif(rmp_kind==3)then
         dp_mom = dp_1(1,1,:)
         dp_s_mom = dp_0(1,1,:)
         ttmp_mom = 0._r8
         ttmp_s_mom = ttmp(1,1,:,1)*dp_0_inv(1,1,:)

         call remapping_core_h(CS,nlev,dp_s_mom,ttmp_s_mom,nlev,dp_mom,ttmp_mom)

         t_1(1,1,:,np1) = ttmp_mom

!         call remapping_core_h(CS,nlev,dp_s_mom,qdp_s_tmp,nlev,dp_mom,qdp_tmp)
         !Tracer qdp
!         elem(ie)%state%Qdp(1,1,:,6,np1_qdp) = qdp_tmp
    else
         call remap1_nofilter(ttmp,np,2,dp_0,dp_1)   
    end if

    if(rmp_kind/=3)then
         t_1(:,:,:,np1)=ttmp(:,:,:,1)*dp_1_inv  !E_rmp
!         q_test_s(:,:,:) = elem(ie)%state%Qdp(:,:,:,6,np1_qdp)/dp_dry(:,:,:)
    end if

    !Remap velocities:

    ttmp(:,:,:,1)=v_0(:,:,1,:,np1)*dp_0 !u*dp_star_moist
    ttmp(:,:,:,2)=v_0(:,:,2,:,np1)*dp_0 !v*dp_star_moist
    if(rmp_kind==1)then
         call remap1(ttmp,np,1,2,2,dp_0,dp_1) !u,v_rmp*dp_moist
    else if(rmp_kind==2)then
         call remap_Q_ppm(ttmp,np,1,2,2,dp_0,dp_1,2)
    elseif(rmp_kind==3)then
         dp_mom = dp_1(1,1,:)
         dp_s_mom = dp_0(1,1,:)
         !remap u
         ttmp_mom = 0._r8
         ttmp_s_mom = ttmp(1,1,:,1)*dp_0_inv(1,1,:)
         call remapping_core_h(CS,nlev,dp_s_mom,ttmp_s_mom,nlev,dp_mom,ttmp_mom)
         ttmp(1,1,:,1) = ttmp_mom
         !remap v
         ttmp_mom = 0._r8
         ttmp_s_mom = ttmp(1,1,:,2)*dp_0_inv(1,1,:)
         call remapping_core_h(CS,nlev,dp_s_mom,ttmp_s_mom,nlev,dp_mom,ttmp_mom)
         ttmp(1,1,:,2) = ttmp_mom
    else
         call remap1_nofilter(ttmp,np,2,dp_0,dp_1)
    end if

    if(rmp_kind/=3)then
            if ( .not. se_prescribed_wind_2d ) &
                 v_1(:,:,1,:,np1)=ttmp(:,:,:,1)*dp_1_inv    !u_rmp
            if ( .not. se_prescribed_wind_2d ) &
                 v_1(:,:,2,:,np1)=ttmp(:,:,:,2)*dp_1_inv    !v_rmp
    endif

    ttmp(:,:,:,1)=t_1(:,:,:,np1) !E_rmp as t
    E_1 = ttmp
 
    !Recover Temperature from Energy (Lin'04):
    if(phi_inc)then
      do k=nlev,1,-1
         num_tmp(:,:,k) = t_1(:,:,k,np1) - &
                          ((v_1(:,:,1,k,np1)**2 + &
                          v_1(:,:,2,k,np1)**2)/2._r8) - phi(:,:,k+1)

         num_tmp(:,:,k) = num_tmp(:,:,k)/cpair
         den_tmp(:,:,k) = 1 - pint_1(k)*kappa*(log(pint_1(k+1))-log(pint_1(k)))!*dp_1_inv(:,:,k)
         t_1(:,:,k,np1) = num_tmp(:,:,k)/den_tmp(:,:,k) !T_rmp              
         phi(:,:,k) = phi(:,:,k+1)+r_universal*t_1(:,:,k,np1)*&
                          (log(pint_1(k+1))-log(pint_1(k))) !phi_rmp
      end do
    else
         t_1(:,:,:,np1) = ( t_1(:,:,:,np1) - (v_1(:,:,1,:,np1)**2 + v_1(:,:,2,:,np1)**2)/2._r8)/cpair  !T_rmp (E_rmp as t is used)
    end if
  
  end subroutine remap_E_cons

  subroutine remap_T_cons(rmp_kind,dp_0,t_0,v_0,pint_0,dp_1,t_1,v_1,pint_1,elem,np1)
    use hybvcoord_mod, only          : hvcoord_t
    use vertremap_mod,          only : remap1, remap1_nofilter,remap_Q_ppm
    use hybrid_mod            , only : hybrid_t!, set_region_num_threads
    use fvm_control_volume_mod, only : fvm_struct
    use control_mod,            only : se_prescribed_wind_2d
    use dimensions_mod        , only : ntrac
    use dimensions_mod        , only : qsize_condensate_loading,qsize_condensate_loading_idx
    use hack_vert_rmp
    use spmd_utils            , only: masterproc
    !use physconst,              only: r_universal
    use MOM_remapping

    type (element_t), intent(inout) :: elem(:)

   !
    type (hvcoord_t) :: hvcoord
    real (kind=r8)   :: dt
    integer          :: ie,i,j,k,nets,nete,np1_qdp,q, m_cnst,temp_int,unitn
    character(len=256):: filename
    real (kind=r8), dimension(np,np,nlev)  :: dp_inv,dp_s_inv
    real (kind=r8), dimension(np,np,nlev)  :: q_test,q_test_s
    real (kind=r8), dimension(np,np,nlev,2):: ttmp
    real (kind=r8), dimension(np,np,nlev+1):: phi
    real(KIND=r8), dimension(nlev+1) :: pint1,pint2,dx1
    real(KIND=r8), dimension(nlev) :: t0,lnp1,lnp2
    logical          :: remap_te,filtered,phi_inc,ppm,pqm
    real(KIND=r8) :: r_universal
    real(KIND=r8) :: dpphi,kappa
    real(KIND=r8),dimension(np,np,nlev) :: num_tmp,den_tmp
    type(remapping_CS)                  :: CS
    real, dimension(nlev) :: dp_mom,ttmp_mom,dp_s_mom,ttmp_s_mom,qdp_tmp,qdp_s_tmp

    integer, intent(in) :: rmp_kind,np1

!    real (kind=r8), dimension(np,np,nlev,2),intent(in) :: E_0
!    real (kind=r8), dimension(np,np,nlev,2),intent(out) :: E_1

    real (kind=r8), dimension(np,np,2,nlev,timelevels),intent(in) :: v_0  !Vel. last arg is 'timelevels'
    real (kind=r8), dimension(np,np,2,nlev,timelevels),intent(inout) :: v_1

    real (kind=r8), dimension(np,np,nlev,timelevels),intent(in) :: t_0  !Temp. last arg same
    real (kind=r8), dimension(np,np,nlev,timelevels),intent(inout) :: t_1

    real (kind=r8), dimension(np,np,nlev),intent(in) :: dp_0,dp_1
    real (kind=r8), dimension(nlev+1),intent(in) :: pint_0,pint_1


!    elem(ie)%state%v = v_0
!    elem(ie)%state%t = t_0
    !ttmp = elem(ie)%state%t(:,:,:,1:2)
!    pint1 = pint_0
!    pint2 = pint_1
!    phi = 0._r8

!    dp_star_moist = dp_0
    dp_s_inv = 1._r8/dp_0
!    dp_moist = dp_1
    dp_inv = 1._r8/dp_1

    
    ttmp(:,:,:,1) = t_0(:,:,:,np1) !T        
!    elem(ie)%state%Qdp(:,:,:,6,np1_qdp) = q_test(:,:,:)*dp_star_dry(:,:,:)

    if(i==1)then
        !    call write_data_TE(pint1,ttmp(1,1,:,1),i)
    end if

    
    ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_0 !T*dp_star_moist

    if(rmp_kind==1)then
            call remap1(ttmp,np,1,1,1,dp_0,dp_1) !T_rmp*dp_moist
    elseif(rmp_kind==2)then
            call remap_Q_ppm(ttmp,np,1,1,1,dp_0,dp_1,1)
!            call remap_Q_ppm(elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp),np,1,qsize,qsize,dp_star_dry,dp_dry,2)
    elseif(rmp_kind==3)then
            dp_mom = dp_1(1,1,:)
            dp_s_mom = dp_0(1,1,:)
            ttmp_mom = 0._r8
            ttmp_s_mom = t_0(1,1,:,np1)

            call remapping_core_h(CS,nlev,dp_s_mom,ttmp_s_mom,nlev,dp_mom,ttmp_mom)

            t_1(1,1,:,np1) = ttmp_mom

            !Remap q only           

!            call remapping_core_h(CS,nlev,dp_s_mom,qdp_s_tmp,nlev,dp_mom,qdp_tmp)
            !Tracer qdp
!            elem(ie)%state%Qdp(1,1,:,6,np1_qdp) = qdp_tmp!*dp_star_dry(1,1,:)

    else
            call remap1_nofilter(ttmp,np,2,dp_0,dp_1)
    end if

    if(rmp_kind/=3)then
            t_1(:,:,:,np1)=ttmp(:,:,:,1)*dp_inv  !T_rmp 
!            q_test_s(:,:,:) = elem(ie)%state%Qdp(:,:,:,6,np1_qdp)/dp_dry(:,:,:)
    end if

        !Here we would print the tracer if we had it working ****


    !Remap velocities:

    ttmp(:,:,:,1)=v_0(:,:,1,:,np1)*dp_0 !u*dp_star_moist
    ttmp(:,:,:,2)=v_0(:,:,2,:,np1)*dp_0 !v*dp_star_moist
    if(rmp_kind==1)then
         call remap1(ttmp,np,1,2,2,dp_0,dp_1) !u,v_rmp*dp_moist
    else if(rmp_kind==2)then
         call remap_Q_ppm(ttmp,np,1,2,2,dp_0,dp_1,2)
    elseif(rmp_kind==3)then
         dp_mom = dp_1(1,1,:)
         dp_s_mom = dp_0(1,1,:)

         !remap u
         ttmp_mom = ttmp(1,1,:,1)
         ttmp_s_mom = ttmp(1,1,:,1)*dp_s_inv(1,1,:)
         call remapping_core_h(CS,nlev,dp_s_mom,ttmp_s_mom,nlev,dp_mom,ttmp_mom)
         v_1(1,1,1,:,np1) = ttmp_mom

         !remap v
         ttmp_mom = ttmp(1,1,:,2)
         ttmp_s_mom = ttmp(1,1,:,2)*dp_s_inv(1,1,:)
         call remapping_core_h(CS,nlev,dp_s_mom,ttmp_s_mom,nlev,dp_mom,ttmp_mom)
         v_1(1,1,2,:,np1) = ttmp_mom
    else
         call remap1_nofilter(ttmp,np,2,dp_0,dp_1)
    end if

    if(rmp_kind/=3)then
            if ( .not. se_prescribed_wind_2d ) &
                 v_1(:,:,1,:,np1)=ttmp(:,:,:,1)*dp_inv    !u_rmp
            if ( .not. se_prescribed_wind_2d ) &
                 v_1(:,:,2,:,np1)=ttmp(:,:,:,2)*dp_inv    !v_rmp
    endif    

    end subroutine remap_T_cons


  subroutine remap_back_forth(num_times,cons_kind,rmp_kind,dp_0,E_0,t_0,v_0,pint_0,dp_1,E_1,t_1,v_1,pint_1,elem,np1) 
  !This routine remaps the same quantity back and forth and generates files for
  !further diagnosis of the remapping algorithm.
  

    use hybvcoord_mod, only          : hvcoord_t
    use vertremap_mod,          only : remap1, remap1_nofilter,remap_Q_ppm
    use hybrid_mod            , only : hybrid_t!, set_region_num_threads
    use fvm_control_volume_mod, only : fvm_struct
    use control_mod,            only : se_prescribed_wind_2d
    use dimensions_mod        , only : ntrac
!    use dimensions_mod        , only : size_condensate_loading,qsize_condensate_loading_idx
    use hack_vert_rmp
    use spmd_utils            , only: masterproc
    !use physconst,              only: r_universal
    use MOM_remapping

    type (element_t), intent(inout) :: elem(:)

   !THESE variables declaration is mostly a copy paste from vertical_remap,
   !needs cleaning
    type (hvcoord_t) :: hvcoord
    real (kind=r8)   :: dt
    integer          :: ie,i,j,k,nets,nete,np1_qdp,q, m_cnst,temp_int,unitn
    character(len=256):: filename
!    real (kind=r8), dimension(np,np,nlev)  :: dp_moist,dp_star_moist,dp_inv,dp_dry,dp_star_dry,dp_s_inv
    real (kind=r8), dimension(np,np,nlev)  :: q_test,q_test_s
    real (kind=r8), dimension(np,np,nlev,2):: ttmp
    real (kind=r8), dimension(np,np,nlev+1):: phi
    real(KIND=r8), dimension(nlev+1) :: pint1,pint2,dx1
    real(KIND=r8), dimension(nlev) :: t0,lnp1,lnp2
    logical          :: remap_te,filtered,phi_inc,ppm,pqm
    real(KIND=r8) :: r_universal
    real(KIND=r8) :: dpphi,kappa
    real(KIND=r8),dimension(np,np,nlev) :: num_tmp,den_tmp
    type(remapping_CS)                  :: CS
    real, dimension(nlev) :: dp_mom,ttmp_mom,dp_s_mom,ttmp_s_mom,qdp_tmp,qdp_s_tmp

    integer, intent(in) :: num_times,cons_kind,rmp_kind,np1

    real (kind=r8), dimension(np,np,nlev,2),intent(inout) :: E_0
    real (kind=r8), dimension(np,np,nlev,2),intent(inout) :: E_1

    real (kind=r8), dimension(np,np,2,nlev,timelevels),intent(inout) :: v_0  !Vel. last arg is 'timelevels'
    real (kind=r8), dimension(np,np,2,nlev,timelevels),intent(inout) :: v_1

    real (kind=r8), dimension(np,np,nlev,timelevels),intent(inout) :: t_0  !Temp. last arg same as v
    real (kind=r8), dimension(np,np,nlev,timelevels),intent(inout) :: t_1

    real (kind=r8), dimension(np,np,nlev),intent(in) :: dp_0,dp_1
    real (kind=r8), dimension(nlev+1),intent(in) :: pint_0,pint_1

    !============
!    cons_kind(s):
    ! 1.- Temperature conserving (no Kinetic or phi energies)
    ! 2.- Energy conserving (according to Lin'04 energy conservation formula)
    !============


    !============
!    rmp_kind(s):
    ! 1.- filtered (remap1_filter())
    ! 2.- ppm (remap_Q_ppm)
    ! 3.- mom (remapping_core_h() using method from vertical_remap() line 1131
    ! -> initialize_remapping())
    ! 4.- no filter (remap1_nofilter())
    !============

    !ttmp = E_0
    !elem(ie)%state%v = v_0
    !elem(ie)%state%t = t_0
    !pint1 = pint_0
    !pint2 = pint_1
    phi = 0._r8
    unitn = 8
    !dp_star_moist = dp_0
    !dp_s_inv = 1._r8/dp_0
    !dp_moist = dp_1
    !dp_inv = 1._r8/dp_1

    

    if(cons_kind==1)then    
        do i = 1,num_times
                call remap_E_cons(rmp_kind,dp_0,E_0,t_0,v_0,pint_0,dp_1,E_1,t_1,v_1,pint_1,elem,np1)
                if(masterproc)then
                        if(i==1)then
                                write (filename, '("TE_eul.dat")' )
                                open(unitn, file=trim(filename),status='replace')
                                write(unitn,*) SUM(t_1(1,1,:,np1)*dp_1(1,1,:)),cpair
                                close(unitn)
                          else
                                write (filename, '("TE_eul.dat")' )
                                open(unitn,file=trim(filename),status='old',position='append')
                                write(unitn,*) SUM(t_1(1,1,:,np1)*dp_1(1,1,:) )
                                close(unitn)
                          end if
                        endif  
        !        if(masterproc)then
        !                if(i==1)then
        !                        write (filename, '("Q_eul_",I0.3,".dat")' ) i
        !                        open(unitn,file=trim(filename),status='replace')
        !                        do k=1,nlev
        !                         if(pqm)then
        !                           write(unitn,*) qdp_tmp(k)
!elem(ie)%state%Qdp(1,1,k,6,np1_qdp)
        !                        else
        !                           write(unitn,*) q_test_s(1,1,k)
        !                         endif
        !                        enddo
        !                        close(unitn)
        !                else
        !                       write (filename, '("Q_eul_",I0.3,".dat")' ) i
        !                       open(unitn,file=trim(filename),status='replace')
        !                       do k=1,nlev
        !                        if(pqm)then
        !                              write(unitn,*) qdp_tmp(k)
!elem(ie)%state%Qdp(1,1,k,6,np1_qdp)
                                               !elem(ie)%state%Qdp(1,1,k,6,np1_qdp)!*dp_moist(1,1,k)
        !                        else
        !                              write(unitn,*) q_test_s(1,1,k)
        !                        endif
        !                       enddo
        !                       close(unitn)
        !        end if
        !                endif
        
                call remap_E_cons(rmp_kind,dp_1,E_1,t_1,v_1,pint_1,dp_0,E_0,t_0,v_0,pint_0,elem,np1)

                !WRITE FILES lagrangian stage i
                call write_data(pint_0,t_0(1,1,:,np1),v_0(1,1,1,:,np1),i)
                call write_data_TE(pint_0,t_0(1,1,:,np1),i)
                call diagnostic(dp_0*v_0(:,:,1,:,np1),&
                                dp_0*t_0(:,:,:,np1)*cpair,&
                                dp_0*((v_0(:,:,1,:,np1)**2 +v_0(:,:,2,:,np1)**2)/2._r8),0._r8,i)

                if(i==num_times)then
                      exit
                endif
        enddo

    elseif(cons_kind==2)then
        do i = 1,num_times
                call remap_T_cons(rmp_kind,dp_0,t_0,v_0,pint_0,dp_1,t_1,v_1,pint_1,elem,np1)

                !WRITE FILES eulerian ?

                call remap_T_cons(rmp_kind,dp_1,t_1,v_1,pint_1,dp_0,t_0,v_0,pint_0,elem,np1)

                !WRITE FILES lagrangian stage i
                call write_data(pint_0,t_0(1,1,:,np1),v_0(1,1,1,:,np1),i)

                call diagnostic(dp_0*v_0(:,:,1,:,np1),&
                                dp_0*t_0(:,:,:,np1)*cpair,&
                                dp_0*((v_0(:,:,1,:,np1)**2 + v_0(:,:,2,:,np1)**2)/2._r8),0._r8,i)


                if(i==num_times)then
                        exit
                endif
        enddo
    
    end if


    !call exit(0) 

   end subroutine remap_back_forth

end module prim_advection_mod
