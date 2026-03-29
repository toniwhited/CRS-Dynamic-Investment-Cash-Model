program evaluate_crs
   !------------------------------------------------------------------------------
   ! PROGRAM: evaluate_crs
   !
   ! PURPOSE:
   !   Solves the CRS model and calculates moments for a given set of parameters.
   !   Simple driver program for testing and comparative statics.
   !
   ! AUTHOR: Adapted from Toni Whited's evaluate_function
   ! DATE: 2025
   !
   !------------------------------------------------------------------------------
   use datatype
   use sizes
   use globals
   use omp_lib
   implicit none

   ! Parameters array
   real(dp) :: params(nop)
   
   ! Moments and value function
   real(dp), allocatable :: momvector(:)
   real(dp), allocatable :: vg(:,:)
   
   ! Timing variables
   integer :: date_time(8), time1(8)
   character(len=12) :: real_clock(3)
   real(dp) :: diter1, diter2
   
   ! Loop variables
   integer :: ii
   
   ! File names
   character(len=60) :: paramfile, momfile
   
   ! Logical for file check
   logical :: c1
   c1 = .true.
   
   !===========================================================================
   ! Start the clock
   !===========================================================================
   call date_and_time(real_clock(1), real_clock(2), real_clock(3), date_time)
   call date_and_time(values=time1)
   diter1 = time1(5)*3600 + time1(6)*60 + time1(7) + time1(8)*0.001_dp
   
   write(*,"('Start date: ',i4,'-',i2.2,'-',i2.2)") date_time(1:3)
   write(*,"('Start time: ',i2.2,':',i2.2,':',i2.2)") date_time(5:7)
   write(*,*)
   
   !===========================================================================
   ! Allocate arrays
   !===========================================================================
   allocate(momvector(nmom))
   allocate(vg(nc, nz))
   allocate(simdata(nyears, nfirms))
   vg = 0.0_dp
   
   !===========================================================================
   ! Set or read parameters
   !===========================================================================
   ! Parameter order:
   !   params(1) = psi      : adjustment cost of investment
   !   params(2) = delta    : depreciation rate
   !   params(3) = mu       : drift of profitability shock
   !   params(4) = rho      : autocorrelation of profitability shock
   !   params(5) = sig_z    : std. dev. of profitability shocks
   
   ! Default parameter values (modify as needed)

   params(1) =  29.9660699147452760_dp    ! psi
   params(2) =   0.0449053817455572_dp    ! delta
   params(3) =  -2.2067113003672532_dp    ! mu
   params(4) =   0.8349248398664271_dp    ! rho
   params(5) =   0.1594425561224784_dp    ! sig_z


   
  ! Optionally read parameters from file
  paramfile = "Input/estfil.txt"
  inquire(file=trim(paramfile), exist=c1)
  if (c1) then
     open(unit=10, file=trim(paramfile), status='old')
     do ii = 1, nop
        read(10, *) params(ii)
     end do
     close(10)
     write(*,*) "Parameters read from file: ", trim(paramfile)
  else
     write(*,*) "Using default parameters"
  end if
   
   ! Print parameters
   write(*,*)
   write(*,*) "====== Parameters ======"
   write(*,"(A20, F10.4)") "psi     = ", params(1)
   write(*,"(A20, F10.4)") "delta   = ", params(2)
   write(*,"(A20, F10.4)") "mu      = ", params(3)
   write(*,"(A20, F10.4)") "rho     = ", params(4)
   write(*,"(A20, F10.4)") "sig_z   = ", params(5)
   write(*,*)
   !===========================================================================
   ! Solve the model
   !===========================================================================
   write(*,*) "====== Solving Model ======"
   write(*,*)
   call shockdraw
   call valfun(params, momvector)
   
   !===========================================================================
   ! Report results
   !===========================================================================
   write(*,*)
   write(*,*) "====== Model Moments ======"
   do ii = 1, nmom
      write(*,"(A10, I2, A3, F12.6)") "Moment ", ii, " = ", momvector(ii)
   end do
   
   ! Write moments to file
   momfile = "Outputcomp/moments.txt"
   open(unit=20, file=trim(momfile), status='replace')
   write(20,*) "====== Parameters ======"
   write(20,"(A20, F10.4)") "psi     = ", params(1)
   write(20,"(A20, F10.4)") "delta   = ", params(2)
   write(20,"(A20, F10.4)") "mu      = ", params(3)
   write(20,"(A20, F10.4)") "rho     = ", params(4)
   write(20,"(A20, F10.4)") "sig_z   = ", params(5)
   write(20,*)
   write(20,*) "====== Moments ======"
   do ii = 1, nmom
      write(20,"(A10, I2, A3, F12.6)") "Moment ", ii, " = ", momvector(ii)
   end do
   close(20)
   
   write(*,*)
   write(*,*) "Results written to: ", trim(momfile)
   
   !===========================================================================
   ! End timing
   !===========================================================================
   call date_and_time(real_clock(1), real_clock(2), real_clock(3), date_time)
   call date_and_time(values=time1)
   diter2 = time1(5)*3600 + time1(6)*60 + time1(7) + time1(8)*0.001_dp
   
   write(*,*)
   write(*,"('End date: ',i4,'-',i2.2,'-',i2.2)") date_time(1:3)
   write(*,"('End time: ',i2.2,':',i2.2,':',i2.2)") date_time(5:7)
   write(*,"('Elapsed time = ',f12.3,' seconds')") diter2 - diter1
   
   !===========================================================================
   ! Deallocate
   !===========================================================================
   deallocate(momvector)
   deallocate(vg)
   
end program evaluate_crs
