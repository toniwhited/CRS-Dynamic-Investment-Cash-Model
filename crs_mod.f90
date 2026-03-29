module datatype
   ! This module defines double precision and makes all of the data structures

   use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64

   implicit none
   integer, parameter :: sp        = real32
   integer, parameter :: dp        = real64

   type :: simulated                     ! simulated data
      sequence
      real(dp) :: v
      real(dp) :: c
      real(dp) :: d
      real(dp) :: i
      real(dp) :: z

   end type simulated

   type :: interp_map
      integer  :: idw   ! lower bracket index
      integer  :: iup   ! upper bracket index
      real(dp) :: wgt   ! weight on upper point: f = wgt*f(iup) + (1-wgt)*f(idw)
   end type interp_map

end module datatype


module sizes
   use datatype
   implicit none

   real(dp), parameter :: pi  = 3.141592653589793238462643383279502884197, maxfunc = 100000.0
   ! Number of grid points
   integer, parameter  :: nz  = 21     !Profitability grid
   integer, parameter  :: nc  = 81     !Cash to capital grid
   integer, parameter  :: ncp = 201    !Cash to capital policy grid
   integer, parameter  :: ni  = 201    ! Number of investment grid values
   integer, parameter :: nhoward = 50  ! Number of Howard policy iterations

   ! Number of parameters and number of moments

   integer, parameter  :: nop = 5 
   integer, parameter  :: allnmom = 7, nmom = 7, numinfl=7 

   ! Specify the parameters in the simulations
   integer, parameter  :: nFirms = 1000
   integer, parameter  :: nYears = 110
   integer, parameter  :: sim_ind = 100

   ! These are the Fortran output file numbers.

   integer, parameter  :: sout=53
   integer, parameter  :: lout=12
   integer, parameter  :: simannout=200

   ! This controls whether you print out the little files

   integer, parameter  :: printfiles=1
   integer, parameter  :: verbose=0


   ! This controls what weight matrix you use and how you calculate standard errors
   !
   !  Complicated = 1 means that you use the (nonoptimal) within weight matrix when you minimize the GMM objective function.
   !                  You then use the clustered weight matrix when calculating the standard errors.
   !
   !  Complicated = 0 means that you just use the clustered weight matrix throughout.
   !
   !
   !  Complicated = 0 puts very little weight on the means. This is bad.
   !
   !  Complicated = 2 means that you use the diagonal within weight matrix when you minimize the GMM objective function.
   !                  You then use the clustered weight matrix when calculating the standard errors.


   integer, parameter  :: complicated = 1

   ! Here we will add some estimation type identifiers that will help
   ! deal with repeating this analysis for early/late and big/small.
   ! The one tricky thing is that the sample size will have to be manually
   ! adjusted for each specification.

   
end module sizes

module globals
   use datatype
   use sizes
   implicit none
   type(simulated), allocatable :: simdata(:,:)                                            ! simulated data
   real(dp), allocatable :: draws(:,:)                                                     ! shock draws for the simulated data
   real(dp), allocatable :: start_i(:), start_c(:)                                         ! Starting investment and cash for the simulation
   integer :: nfcnev=0
   integer :: errcode ! This is the number of simulated anneal evalutations and an error code
   real(dp) :: lambda = 0.0542355_dp ! From GWZ (2021) RFS
end module globals


module pickmoments
   use datatype
   use sizes

   implicit none




   !=======================note: this is an abbreviated moment list==============================================================
   integer, parameter,dimension(nmom)  :: pickout = (/ 1, 2, 3, 4, 5, 6, 7 /)

   integer, parameter,dimension(numinfl) :: isvar = 0


   !===================Make cryptic moment names.================================
   character(len=11),parameter,dimension(allnmom) :: amomname = (/     &
   !
   !
   !
   !                                       ! This is a list of the possible moments.
   !
      'mu___cash  ',   &       !allmoms(1)
      'v____cash  ',   &       !allmoms(2)
      'mu___invest',   &       !allmoms(3)
      'v____invest',   &       !allmoms(4)
      'mu___profit',   &       !allmoms(5)
      'v____profit',   &       !allmoms(6)
      'rho__profit'    /)      !allmoms(7)





   character (len = 11), parameter, dimension(nop) :: pname = (/    &
      'psi        ',   &      ! pname(1) =
      'delta      ',   &      ! pname(2) =
      'mu         ',   &      ! pname(3) =
      'rho        ',   &      ! pname(4) =
      'sigma      '   /)      ! pname(5) =

end module pickmoments

module myseed
   use datatype
   implicit none
   integer, allocatable :: iiseed(:)
   integer :: iiiz, ifi, iyi, nnn


contains
   subroutine init_seed()
      intrinsic random_seed, random_number
      call random_seed(size = nnn)
      allocate(iiseed(nnn))
      iiseed(1) = 987654321
      do iiiz=2,nnn
         iiseed(iiiz) = iiseed(iiiz-1)+1
      enddo
      call random_seed(PUT=iiseed)
   end subroutine
end module myseed
