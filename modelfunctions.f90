subroutine makegrids(params, ig, cg, cpg)
   !------------------------------------------------------------------------------
   ! SUBROUTINE: makegrids
   !
   ! PURPOSE:
   !   Constructs the investment and cash/debt grids for the CRS model.
   !
   ! ARGUMENTS:
   !   delta    [IN]  : Depreciation rate
   !   theta    [IN]  : Collateral parameter (determines cash grid bounds)
   !   delp     [IN]  : Grid spacing parameter for investment
   !   ig(:)    [OUT] : Investment grid, dimension (ni)
   !   cg(:)    [OUT] : Cash/debt state grid, dimension (nc)
   !   cpg(:)   [OUT] : Cash/debt policy grid, dimension (ncp)
   !
   ! GRID SPECIFICATIONS:
   !   Investment: Centered around delta, with spacing delta/delp
   !   Cash/Debt:  Uniform grid from 0 to theta
   !
   ! AUTHOR: Adapted from Toni Whited's CRS model
   ! DATE: 2025
   !
   !------------------------------------------------------------------------------
   use datatype
   use sizes
   implicit none

   ! Arguments
   real(dp), intent(in)  :: params(nop)
   real(dp), intent(out) :: ig(ni)
   real(dp), intent(out) :: cg(nc)
   real(dp), intent(out) :: cpg(ncp)

   ! Local variables
   real(dp) :: mini, max_c, min_c, delp = 4.0_dp
   real(dp) :: delta, theta
   integer  :: ii, ic

   theta = 0.7_dp ! upper bound for the cash grid. 
   delta = params(2)

   !==========================================================================
   ! Investment Grid
   !==========================================================================
   ! Grid is centered around delta with spacing delta/delp
   ! mini is the minimum investment rate

   mini = delta * (1.0_dp - real(ceiling(real(ni) / (2.0_dp * delp)))) + delta

   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
   !$OMP DO SCHEDULE(dynamic)
   do ii = 1, ni
      ig(ii) = mini + real(ii - 1) * delta / delp
   end do
   !$OMP END DO
   !$OMP END PARALLEL

   !==========================================================================
   ! Cash/Debt State Grid
   !==========================================================================
   ! Uniform grid from -theta to +theta
   max_c = theta
   min_c = 0.0_dp

   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ic)
   !$OMP DO SCHEDULE(dynamic)
   do ic = 1, nc
      cg(ic) = min_c + real(ic - 1) * (max_c - min_c) / real(nc - 1)
   end do
   !$OMP END DO
   !$OMP END PARALLEL
 
   !==========================================================================
   ! Cash/Debt Policy Grid
   !==========================================================================
   ! Uniform grid from -theta to +theta (can have different resolution than cg)
   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ic)
   !$OMP DO SCHEDULE(dynamic)
   do ic = 1, ncp
      cpg(ic) = min_c + real(ic - 1) * (max_c - min_c) / real(ncp - 1)
   end do
   !$OMP END DO
   !$OMP END PARALLEL

end subroutine makegrids

subroutine printstatespaces(cg, ig, cpg, zg)
   use datatype
   use sizes
   implicit none

   real(dp), intent(in) :: cg(nc)    ! debt grid
   real(dp), intent(in) :: ig(ni)    ! investment grid
   real(dp), intent(in) :: cpg(ncp)  ! debt policy grid
   real(dp), intent(in) :: zg(nz)    ! shock grid

   integer :: iunit

   iunit = 77
   open(iunit, file='Outputcomp/statespaces.txt', status='replace')

   write(iunit,*) '============================================='
   write(iunit,*) '           State Space Grids'
   write(iunit,*) '============================================='
   write(iunit,*) ''

   write(iunit,*) 'Cash grid (cg), size =', nc
   write(iunit,'(f12.4)') cg

   write(iunit,*) ''
   write(iunit,*) 'Investment grid (ig), size =', ni
   write(iunit,'(f12.4)') ig

   write(iunit,*) ''
   write(iunit,*) 'Cash policy grid (cpg), size =', ncp
   write(iunit,'(f12.4)') cpg

   write(iunit,*) ''
   write(iunit,*) 'Shock grid (zg), size =', nz
   write(iunit,'(f12.4)') zg

   write(iunit,*) ''
   write(iunit,*) '============================================='

   close(iunit)

   open(78, file='Outputcomp/cg.txt', status='replace')
   write(78,'(f12.4)') cg
   close(78)
   open(79, file='Outputcomp/zg.txt', status='replace')
   write(79,'(f12.4)') zg
   close(79)

end subroutine printstatespaces

subroutine bellman_maximize(beeta, payoff, queue_i, vnew, gidx_c, gidx_i)
   !
   ! Performs the Bellman maximization step
   !
   use datatype
   use sizes
   use omp_lib
   implicit none

   real(dp), intent(in) :: beeta
   real(dp), intent(in) :: payoff(nc,ncp,ni,nz)
   real(dp), intent(in) :: queue_i(ncp,ni,nz)
   real(dp), intent(inout) :: vnew(nc,nz)
   integer, intent(inout) :: gidx_c(nc,nz)
   integer, intent(inout) :: gidx_i(nc,nz)

   integer :: ic, iic, iik, iz
   real(dp) :: xtemp

   vnew = -huge(0.0_dp)
   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(xtemp,iic,iik,ic,iz)
   do iic=1,ncp
      do iik=1,ni
         !$OMP DO COLLAPSE(2)
         do ic=1,nc
            do iz=1,nz
               xtemp = beeta*queue_i(iic,iik,iz)+payoff(ic,iic,iik,iz)
               if (xtemp > vnew(ic,iz)) then
                  vnew(ic,iz) = xtemp
                  gidx_c(ic,iz) = iic
                  gidx_i(ic,iz) = iik
               endif
            enddo
         enddo
         !$OMP END DO
      enddo
   enddo
   !$OMP END PARALLEL

end subroutine bellman_maximize


subroutine howard_improve(beeta, payoff, queue_i, vnew, gidx_c, gidx_i)
   !
   ! Performs Howard policy improvement iterations:
   ! Given fixed policy indices, iterates the value function
   ! to accelerate convergence
   !
   use datatype
   use sizes
   use omp_lib
   implicit none

   real(dp), intent(in) :: beeta
   real(dp), intent(in) :: payoff(nc,ncp,ni,nz)
   real(dp), intent(in) :: queue_i(ncp,ni,nz)
   real(dp), intent(inout) :: vnew(nc,nz)
   integer, intent(inout) :: gidx_c(nc,nz)
   integer, intent(inout) :: gidx_i(nc,nz)

   integer :: ic, iz


   ! Interpolate to policy grid if grids differ

   ! Evaluate at fixed policy (no maximization)
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ic,iz) COLLAPSE(2) SCHEDULE(static)
   do iz = 1,nz
      do ic = 1,nc
         vnew(ic,iz) = beeta*queue_i(gidx_c(ic,iz),gidx_i(ic,iz),iz) + payoff(ic,gidx_c(ic,iz),gidx_i(ic,iz),iz)
      enddo
   enddo
   !$OMP END PARALLEL DO

end subroutine howard_improve


subroutine build_interp_map(src_grid, nsrc, tgt_grid, ntgt, imap)
   !
   ! Builds an interpolation map from src_grid to tgt_grid.
   ! For each target point, finds the bracketing source interval
   ! and computes the weight on the upper point.
   ! Works for any grid (uniform or not).
   !
   use datatype
   implicit none

   integer,  intent(in)  :: nsrc, ntgt
   real(dp), intent(in)  :: src_grid(nsrc), tgt_grid(ntgt)
   type(interp_map), intent(out) :: imap(ntgt)

   integer  :: ic, ii, idw
   real(dp) :: x

   do ic = 1, ntgt
      x = tgt_grid(ic)

      ! Search from top for bracketing interval
      idw = 1
      do ii = nsrc, 1, -1
         if (x >= src_grid(ii)) then
            idw = ii
            exit
         endif
      enddo
      idw = min(idw, nsrc - 1)

      imap(ic)%idw = idw
      imap(ic)%iup = idw + 1
      imap(ic)%wgt = (x - src_grid(idw)) / (src_grid(idw + 1) - src_grid(idw))
   enddo

end subroutine build_interp_map


subroutine bellman_setup(v, tmat, delta, ig, cmap, queue_i)
   !
   ! Performs the Bellman setup step:
   ! - Computes expectation: queue = v * tmat
   ! - Interpolates from state grid to policy grid using precomputed map
   ! - Scales by (1 - delta + investment) for each investment level
   !
   use datatype
   use sizes
   use omp_lib
   implicit none

   real(dp), intent(in) :: v(nc,nz)
   real(dp), intent(in) :: tmat(nz,nz)
   real(dp), intent(in) :: delta
   real(dp), intent(in) :: ig(ni)
   type(interp_map), intent(in) :: cmap(ncp)
   real(dp), intent(out) :: queue_i(ncp,ni,nz)

   integer :: iz, ic, ii, idw, iup
   real(dp) :: wgt
   real(dp) :: queue(nc,nz), queuelong(ncp,nz)

   ! Compute expectation
   !$OMP WORKSHARE
   queue = matmul(v,tmat)
   !$OMP END WORKSHARE

   ! Interpolate to policy grid using precomputed map
   if (nc == ncp) then
      queuelong(1:nc,:) = queue
   else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,ic,idw,iup,wgt) COLLAPSE(2) SCHEDULE(static)
      do iz = 1,nz
         do ic = 1,ncp
            idw = cmap(ic)%idw
            iup = cmap(ic)%iup
            wgt = cmap(ic)%wgt
            queuelong(ic,iz) = (1.0_dp - wgt)*queue(idw,iz) + wgt*queue(iup,iz)
         end do
      end do
      !$OMP END PARALLEL DO
   endif

   ! Scale by (1 - delta + investment) for each investment level
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(static)
   do ii = 1,ni
      queue_i(:,ii,:) = queuelong*(1.0_dp - delta + ig(ii))
   end do
   !$OMP END PARALLEL DO

end subroutine bellman_setup

subroutine printpolicies(v, c, i, d, cg, zg)
   !------------------------------------------------------------------------------
   ! SUBROUTINE: printpolicies
   !
   ! PURPOSE:
   !   Writes the value function and policy functions to output files.
   !
   ! ARGUMENTS:
   !   v(nc,nz)   [IN] : Value function
   !   c(nc,nz)   [IN] : Optimal debt policy
   !   i(nc,nz)   [IN] : Optimal investment policy
   !   d(nc,nz)   [IN] : Optimal dividend policy
   !   cg(nc)     [IN] : Cash/debt state grid
   !   ig(ni)     [IN] : Investment grid
   !   cpg(ncp)   [IN] : Cash/debt policy grid
   !   zg(nz)     [IN] : Shock grid
   !
   !------------------------------------------------------------------------------
   use datatype
   use sizes
   implicit none

   real(dp), intent(in) :: v(nc,nz)
   real(dp), intent(in) :: c(nc,nz)
   real(dp), intent(in) :: i(nc,nz)
   real(dp), intent(in) :: d(nc,nz)
   real(dp), intent(in) :: cg(nc)
   real(dp), intent(in) :: zg(nz)

   integer :: iz, ic, iunit

      iunit = 78
      open(iunit, file='Outputcomp/policies.txt', status='replace')

      write(iunit,*) '============================================='
      write(iunit,*) '           Model Answer'
      write(iunit,*) '============================================='
      write(iunit,*) ''

      write(iunit,*) 'Grid dimensions:'
      write(iunit,*) '  nc  (cash state grid)  =', nc
      write(iunit,*) '  ncp (cash policy grid) =', ncp
      write(iunit,*) '  ni  (investment grid)  =', ni
      write(iunit,*) '  nz  (shock grid)       =', nz
      write(iunit,*) ''

      write(iunit,*) '============================================='
      write(iunit,*) '           Value Function v(c,z)'
      write(iunit,*) '============================================='
      write(iunit,*) 'Rows = cash states, Columns = shock states'
      write(iunit,'(a12,100f12.4)') 'z ->', zg
      do ic = 1, nc
         write(iunit,'(f12.3,100f12.4)') cg(ic), v(ic,:)
      end do
      write(iunit,*) ''

      write(iunit,*) '============================================='
      write(iunit,*) '           Cash Policy c''(c,z)'
      write(iunit,*) '============================================='
      write(iunit,*) 'Rows = cash states, Columns = shock states'
      write(iunit,'(a12,100f12.4)') 'z ->', zg
      do ic = 1, nc
         write(iunit,'(f12.3,100f12.4)') cg(ic), c(ic,:)
      end do
      write(iunit,*) ''

      write(iunit,*) '============================================='
      write(iunit,*) '           Investment Policy i(c,z)'
      write(iunit,*) '============================================='
      write(iunit,*) 'Rows = cash states, Columns = shock states'
      write(iunit,'(a12,100f12.4)') 'z ->', zg
      do ic = 1, nc
         write(iunit,'(f12.3,100f12.4)') cg(ic), i(ic,:)
      end do
      write(iunit,*) ''

      write(iunit,*) '============================================='
      write(iunit,*) '           Dividend Policy d(c,z)'
      write(iunit,*) '============================================='
      write(iunit,*) 'Rows = cash states, Columns = shock states'
      write(iunit,'(a12,100f12.4)') 'z ->', zg
      do ic = 1, nc
         write(iunit,'(f12.3,100f12.4)') cg(ic), d(ic,:)
      end do
      write(iunit,*) ''

      write(iunit,*) '============================================='
      close(iunit)

      ! Also write separate files for easy import
      open(4, file='Outputcomp/v.txt', status='replace')
      open(3, file='Outputcomp/c.txt', status='replace')
      open(7, file='Outputcomp/d.txt', status='replace')
      open(8, file='Outputcomp/i.txt', status='replace')

      do iz = 1, nz
         write(4, '(f10.4,6x,100f16.3)') zg(iz), v(:,iz)
         write(3, '(f10.4,6x,100f16.3)') zg(iz), c(:,iz)
         write(7, '(f10.4,6x,100f16.3)') zg(iz), d(:,iz)
         write(8, '(f10.4,6x,100f16.3)') zg(iz), i(:,iz)
      end do
      close(4)
      close(3)
      close(7)
      close(8)

end subroutine printpolicies

subroutine tauchen(rho_z,sigma2_z,mu,ennz,m,z,trans)
   !------------------------------------------------------------------------------
   ! DESCRIPTION:
   !   This subroutine implements the Tauchen method to compute the Markov
   !   transition matrix for an autoregressive process.
   !
   ! MODEL SPECIFICATION:
   !   z - mu = ρ(z(-1) - mu) + ε
   !   where:
   !     - mu is the unconditional mean
   !     - ε ~ N(0, sigma2_z)
   !
   !
   ! AUTHOR: Toni Whited
   ! DATE: 2004
   ! MODIFIED: 2022
   !
   !------------------------------------------------------------------------------
  use datatype
  implicit none
  real(dp), intent(in)  :: rho_z, sigma2_z, m, mu   ! serial correlation and variances of the shock, # std. dev. 
  integer, intent(in)   :: ennz                     ! dimension of shock
  real(dp), intent(out) :: trans(ennz,ennz), z(ennz)! transition matrix and grids
  real(dp) :: sig_z, sig_eps                         ! unconditional and innovation std devs
  real(dp) :: w, binhi, binlo                       ! width of the Tauchen grid, bin endpoints
  integer :: iz, izz
  ! z grid (built centered at 0; shifted by mu after transitions are computed)
  sig_eps = sqrt(sigma2_z)
  sig_z = sqrt(sigma2_z /(1.0d0 - rho_z**2.0d0))
  z = 0.0d0
  do iz = 1, ennz
      z(iz) = -m * sig_z + (2.0d0 * m * sig_z * dble(iz - 1)) / dble(ennz - 1)
  enddo
  trans = 0.0d0
  w = (z(2) - z(1))
  do iz = 1, ennz
      do izz = 1, ennz
          binhi = (z(izz) - rho_z * z(iz) + w / 2.0d0) / sig_eps
          binlo = (z(izz) - rho_z * z(iz) - w / 2.0d0) / sig_eps
          if (izz .eq. 1) then
              trans(iz,izz) = 0.5d0 * (1.0d0 + derf(binhi / sqrt(2.0d0)))
          else if (izz .eq. ennz) then
              trans(iz,izz) = 1.0d0 - (0.5d0 * (1.0d0 + derf(binlo / sqrt(2.0d0))))
          else
              trans(iz,izz) = 0.5d0 * (1.0d0 + derf(binhi / sqrt(2.0d0))) - (0.5d0 * (1.0d0 + derf(binlo / sqrt(2.0d0))))
          endif
      enddo
  enddo
  z = z + mu
end subroutine tauchen