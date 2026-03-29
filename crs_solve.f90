subroutine valfun(params,momvector)
   use datatype
   use sizes
   use globals
   use omp_lib
   implicit none

   real(dp), intent(in) :: params(nop)
   real(dp), intent(out) :: momvector(nmom)



   real(dp), parameter :: toler = 0.0001_dp
   integer, parameter :: capt = 3000

   real(dp) :: mu, sig_z, rho, rf, delta, psi ! parameters

   real(dp) :: str=2.5_dp, sig2z                                   ! stuff to define the grids

   real(dp) :: diff = 100.0_dp, tau_c = 0.10_dp, beeta

   integer :: ii, n, ic, icp, iz, counter, iih
   integer :: pdiff                                                                                  ! Stuff for policy function convergence
   real(dp), allocatable :: tmat(:,:), zg(:)                                                         ! stuff for the transition matrix
   real(dp), allocatable :: ig(:), cg(:), cpg(:)                                                     ! state variable and choice variable grids
   real(dp), allocatable :: payoff(:,:,:,:), vnew(:,:), v(:,:), &                    ! profits, dividends, values
      queue_i(:,:,:), vg(:,:)

   integer, allocatable :: gidx_i(:,:), gidx_c(:,:), oidx_i(:,:), oidx_c(:,:)                                                  ! This is the policy index function!
   real(dp), allocatable :: i(:,:),c(:,:),d(:,:)
   type(interp_map), allocatable :: cmap(:)                                                            ! interpolation map from cg to cpg





   !=========================Allocate all of the arrays=================================
   allocate(zg(nz))
   allocate(tmat(nz,nz))
   allocate(cg(nc))
   allocate(ig(ni))
   allocate(cpg(ncp))
   allocate(payoff(nc,ncp,ni,nz))
   allocate(v(nc,nz))
   allocate(vnew(nc,nz))
   allocate(vg(nc,nz))
   allocate(queue_i(ncp,ni,nz))
   allocate(gidx_i(nc,nz))
   allocate(gidx_c(nc,nz))
   allocate(oidx_i(nc,nz))
   allocate(oidx_c(nc,nz))
   allocate(c(nc,nz))
   allocate(i(nc,nz))
   allocate(d(nc,nz))
   allocate(cmap(ncp))

   !===========================================================================================================
   !===========================================================================================================
   !===========================================================================================================
   !         This section reads in the parameters and sets up the state spaces
   !===========================================================================================================
   !===========================================================================================================
   !===========================================================================================================


   rf       = 0.02_dp
   psi      = params(1)       !adjustment cost of investment.
   delta    = params(2)       !depreciation rate
   mu       = params(3)       !drift of profitability shock
   rho      = params(4)       !autocorrelation of profiability shock
   sig_z    = params(5)       !std. dev. of profitability shocks
   beeta = 1.0_dp/(1.0_dp+rf*(1.0_dp+tau_c))! This approximates the tax benefit of debt

   !===========================================================================================================
   !       This part sets up the shock grid and the transition matrix
   !===========================================================================================================

   n = nz
   sig2z = sig_z**2.0_dp
   call  tauchen(rho,sig2z,mu,n,str,zg,tmat)

   !do iz=1,nz
   !   write (*,"(100f8.4)") tmat(iz,:)
   !enddo

   tmat = transpose(tmat)

   zg = exp(zg)

   !===========================================================================================================
   !       This part sets up the investment and debt grids
   !===========================================================================================================

   call makegrids(params, ig, cg, cpg)
   call build_interp_map(cg, nc, cpg, ncp, cmap)
   
   if (verbose==1) write(*,*) 'Finished making grids.'
   if (printfiles==1)call printstatespaces(cg, ig, cpg, zg)
   !===========================================================================================
   ! Now make the profit flow and initialize everything for the value function iteration


   ! Define the sources and uses of funds identity
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii,icp,iz,ic) COLLAPSE(4) SCHEDULE(static)
   do iz = 1,nz
      do ii = 1,ni
         do icp = 1,ncp
            do ic = 1,nc
               payoff(ic,icp,ii,iz) = zg(iz) - ig(ii) - 0.5_dp*psi*(ig(ii)**2.0_dp) + cg(ic)*(1.0_dp + rf) - cpg(icp)*(1.0_dp - delta + ig(ii))
            enddo
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   if (verbose==1) write(*,*) 'Finished making profit.'
   where (payoff < 0.0_dp)
      payoff = payoff*(1.0_dp+lambda)
   end where

   vg = 0.0_dp
   v = vg
   gidx_c = 1
   gidx_i = 1
   oidx_c = 1
   oidx_i = 1
   errcode = 0
   diff = 100.0_dp
   pdiff = 100
   counter = 0

   !======================================value function iteration=============================
   convergeloop: do

      ! Save old policy indices
      oidx_c = gidx_c
      oidx_i = gidx_i

      ! Bellman setup: expectation, interpolation, scaling
      call bellman_setup(v, tmat, delta, ig, cmap, queue_i)
      call bellman_maximize(beeta, payoff, queue_i, vnew, gidx_c, gidx_i)

      do iih=1,nhoward
         call bellman_setup(vnew, tmat, delta, ig, cmap, queue_i)
         call howard_improve(beeta, payoff, queue_i, vnew, gidx_c, gidx_i)
      enddo

      diff = maxval(abs((v-vnew)))
      pdiff = sum(abs(oidx_i-gidx_i))+sum(abs(oidx_c-gidx_c))
      if (verbose==1)write (*,"('Iteration: ',i6,'  Value Function Diff: 'f12.9,'  Policy Index Diff: ',i6)") counter, diff, pdiff

      !$OMP WORKSHARE
      v=vnew
      !$OMP END WORKSHARE
      counter = counter + 1

      if (diff > 1.0e16)  errcode = 2
      if (counter > capt) errcode = 5

      if ((diff < toler .and. pdiff < 1) .or. errcode > 0 ) then
         exit convergeloop
      endif

   enddo convergeloop
   if (mod(real(nfcnev),100.)==0.) then
      write(*,"('Error Code = ',i4,3x,'Iterations = ',i4,3x,'SA iteration = ',i8)") errcode, counter, nfcnev
   endif


   if (errcode == 0) then
      ! Construct the policy functions from the policy index functions.

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iz,ic)
      !$OMP DO SCHEDULE(dynamic)
      do iz=1,nz
         do ic=1,nc

            c(ic,iz) =  cpg(gidx_c(ic,iz))
            i(ic,iz) =  ig(gidx_i(ic,iz))
            d(ic,iz) =  payoff(ic,gidx_c(ic,iz),gidx_i(ic,iz),iz)

         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      if (printfiles==1) call printpolicies(v, c, i, d, cg, zg)

      momvector = 0.0_dp
      
      call simpanel(v,c,i,d,cg,zg,tmat)

      call makemoments(momvector)
      !=========note that you only update the guess if the model solves!!!!!
      !$OMP WORKSHARE
      vg = vnew
      !$OMP END WORKSHARE
   else
      momvector = -10.0_dp
   endif




   !=========================Deallocate all of the arrays=================================
   deallocate(zg)
   deallocate(tmat)
   deallocate(cg)
   deallocate(ig)
   deallocate(cpg)
   deallocate(payoff)
   deallocate(v)
   deallocate(vnew)
   deallocate(vg)
   deallocate(queue_i)
   deallocate(gidx_i)
   deallocate(gidx_c)
   deallocate(oidx_i)
   deallocate(oidx_c)
   deallocate(c)
   deallocate(i)
   deallocate(d)
   deallocate(cmap)

end subroutine valfun