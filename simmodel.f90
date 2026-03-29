subroutine simpanel(v,c,i,d,cg,zg,tmat)
!==============================================================================
! SUBROUTINE: simpanel
! PURPOSE: Simulates a panel of firms using policy functions and transition
!          matrices to generate artificial firm-level dat
!==============================================================================
      use datatype
      use sizes
      use globals
 
      implicit none

      real(dp), intent(in) :: tmat(nz,nz)
      real(dp), intent(in) :: v(nc,nz)   
      real(dp), intent(in) :: c(nc,nz)
      real(dp), intent(in) :: i(nc,nz)
      real(dp), intent(in) :: d(nc,nz)
      real(dp), intent(in) :: cg(nc)
      real(dp), intent(in) :: zg(nz)


      real(dp), allocatable :: cdf_wgt(:), phatcdf(:,:)
      real(dp), allocatable :: cstart(:), istart(:)
      integer,  allocatable :: ls(:,:), zls(:,:) 

      real(dp) :: erg_dist(nz,nz)
      real(dp) :: vold, iold, cold, zold, dold, vprime, iprime, cprime, dprime, zprime, diffi 
      
      integer :: iz, ifi, iti, picki, pickc, pickz 
 

      !-----------------------------------Now simulate the model.-------------------------
      !
      allocate (phatcdf(nz,nz))
      allocate (cdf_wgt(nz))
      allocate (ls(nyears+1,nfirms))
      allocate (zls(nyears+1,nfirms))
      allocate (istart(nfirms))
      allocate (cstart(nfirms))

      ! Unconditional cdf of the shocks
      ! Pull it from the ergodic distribution

      erg_dist = transpose(tmat)
      do iti = 1,100
         erg_dist = matmul(erg_dist,transpose(tmat))
      enddo


      cdf_wgt = erg_dist(1,:)! 1.0_dp/dble(nz)
      do iz = 2,nz
          cdf_wgt(iz) = cdf_wgt(iz) + cdf_wgt(iz-1)
      enddo

      ! Conditional shock cdf
      phatcdf = transpose(tmat)
      do iz=2,nz
        phatcdf(:,iz) = phatcdf(:,iz) + phatcdf(:,iz-1) 
      enddo

      ! write out neat rows of phatcdf for checking
      !open (unit=261,file='Outputcomp/phatcdf.txt')

      ! This just randomly assigns the starting locations from the unconditional distribution
      do ifi = 1,nfirms
          ls(1,ifi) = nz
          oneloop: do iz = 1,nz
            diffi = draws(1,ifi) - cdf_wgt(iz);
            if (diffi < 0.0_dp) then
               ls(1,ifi) = iz
               exit oneloop
            endif
          enddo oneloop
      enddo

      do ifi = 1,nfirms
         cstart(ifi) = start_c(ifi)
         istart(ifi) = start_i(ifi)
      enddo

      crosssection: do ifi=1,nfirms
        


         picki = min(int(ni*istart(ifi))+1, ni)
         pickc = min(int(nc*cstart(ifi))+1, nc) 
         pickz = ls(1,ifi)

         vold = v(pickc,pickz)
         iold = i(pickc,pickz)
         cold = c(pickc,pickz)
         dold = d(pickc,pickz)

         zls(1,ifi) = ls(1,ifi)

         zold = zg(zls(1,ifi))


         timeseries: do iti=1,nyears
            ls(iti+1,ifi) = nz   !Initialize it at the end in case the loop never gets there

            pickloop: do iz=1,nz
              diffi = draws(iti+1, ifi) - phatcdf(ls(iti, ifi),iz)
              if (diffi < 0.0_dp) then
                 ls(iti+1,ifi) = iz
                 exit pickloop
              endif
            enddo pickloop

            call interpol(ls(iti,ifi),cold,cg,v,vprime)
            call interpol(ls(iti,ifi),cold,cg,i,iprime)
            call interpol(ls(iti,ifi),cold,cg,c,cprime)
            call interpol(ls(iti,ifi),cold,cg,d,dprime)

            zls(iti+1,ifi) = ls(iti+1,ifi)

            zprime = zg(zls(iti+1,ifi))

            simdata(iti,ifi)%v = vprime
            simdata(iti,ifi)%i = iprime
            simdata(iti,ifi)%c = cprime
            simdata(iti,ifi)%z = zprime
            simdata(iti,ifi)%d = dprime

            vold=vprime
            iold=iprime
            cold=cprime
            zold=zprime
            dold=dprime

         enddo timeseries
      enddo  crosssection
       
      if (printfiles == 1) then
            open (unit=262,file='Outputcomp/simdata.txt')
 
         do ifi=1,nfirms
            do iti=1,nyears
                 write(262,"(2i6,20f18.3)") ifi,iti,simdata(iti,ifi)%v,simdata(iti,ifi)%c,simdata(iti,ifi)%i,simdata(iti,ifi)%d,simdata(iti,ifi)%z 
            enddo
         enddo
         close(262)
      endif

      deallocate(phatcdf)
      deallocate(cdf_wgt)
      deallocate(ls)
      deallocate(zls)
      deallocate(istart)
      deallocate(cstart)

end subroutine simpanel

subroutine interpol(idz,cold,cg,v,vprime)


   use datatype
   use sizes
   implicit none
   integer,  intent(in) :: idz
   real(dp), intent(in) :: cold, cg(nc), v(nc,nz)
   real(dp), intent(out) :: vprime
   real(dp) ::  cfrac 
   integer :: idchi, idclo, enn

   enn = size(cg);  
   call pick(cold,cg,idclo,idchi,enn,cfrac)

   !==========================================================================

   vprime = cfrac*v(idchi,idz) + (1.0_dp - cfrac)*v(idclo,idz)  


end subroutine interpol


subroutine pick(iold,kg,idklo,idkhi,enn,kfrac)
!==============================================================================
   use sizes
   use datatype
   use globals
   implicit none
   integer, intent(in)   :: enn
   real(dp), intent(in)  :: iold, kg(enn)
   real(dp), intent(out) :: kfrac
   integer,  intent(out) :: idklo, idkhi
   real(dp) ::  klo, khi
   real(dp), parameter :: eps = 1.0e-6
   integer :: ii

   idklo = 1
   do ii=enn,1,-1
     if (abs(iold-kg(ii))<eps) then
        idklo = ii
        exit
     endif
     if (iold > kg(ii)) then
        idklo = ii
        exit
     endif
   enddo
   if (idklo < enn) then
     idkhi = idklo + 1
     klo = kg(idklo)
     khi = kg(idkhi)
     kfrac = abs(iold-klo)/abs(khi-klo)
   else
      idkhi = enn
      kfrac = 1.0_dp
   endif
end subroutine pick