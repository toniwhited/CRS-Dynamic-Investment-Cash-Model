subroutine shockdraw
!------------------------------------------------------------------------------
! SUBROUTINE: shockdraw
!
! PURPOSE:
!   Generates random shock values for each year and firm, initializing arrays
!   for shocks (`draws`) and starting values for capital (`start_k`) and bonds (`start_b`).
!
! AUTHOR: Toni Whited
! DATE: 2004
! UDAPTED: 2024
!
! ARGUMENTS:
!   None.
!   
! NOTE:
!   Relies on Fortran’s `random_seed` and `random_number` to generate reproducible shocks.
!------------------------------------------------------------------------------

   use datatype
   use sizes
   use globals
   use myseed
   implicit none

   real(dp) ::  shocktemp
   integer :: iyear, ifirm
    call init_seed()   

   allocate (draws(nYears+1,nfirms))
   allocate (start_i(nfirms))
   allocate (start_c(nfirms))

   do iyear = 1,nYears+1
      do ifirm = 1,nfirms
         call random_number(harvest=shocktemp)
         draws(iyear,ifirm) = shocktemp
      enddo
   enddo

   do ifirm = 1,nfirms
      call random_number(harvest=shocktemp)
      start_i(ifirm) = shocktemp
      call random_number(harvest=shocktemp)
      start_c(ifirm) = shocktemp
   enddo
end subroutine shockdraw
