subroutine makemoments(momvector)
    use datatype
    use sizes
    use globals
    use pickmoments
    implicit none
    real(dp), intent(inout) :: momvector(nmom)

    real(dp), allocatable ::  vs(:,:), cs(:,:), ds(:,:), is(:,:), pis(:,:), zs(:,:)
    real(dp), allocatable ::  vms(:,:), cms(:,:), dms(:,:), ims(:,:), pims(:,:), pifs(:,:), cifs(:,:)



    integer :: nobs
    real(dp) :: mu_z, var_z, rho_z
    real(dp) :: mu_i, var_i
    real(dp) :: mu_c, var_c, rho_c 
    real(dp) :: allmom(numinfl)

    allocate(vs(nYears,nFirms))
    allocate(cs(nYears,nFirms))
    allocate(ds(nYears,nFirms))
    allocate(is(nYears,nFirms))
    allocate(pis(nYears,nFirms))
    allocate(zs(nYears,nFirms))

    allocate(vms(nYears-sim_ind+1,nFirms))
    allocate(cms(nYears-sim_ind+1,nFirms))
    allocate(dms(nYears-sim_ind+1,nFirms))
    allocate(ims(nYears-sim_ind+1,nFirms))
    allocate(pims(nYears-sim_ind+1,nFirms))
    allocate(pifs(nYears-sim_ind+1,nFirms))
    allocate(cifs(nYears-sim_ind+1,nFirms))

    vs     =  simdata%v
    cs     =  simdata%c
    ds     =  simdata%d
    is     =  simdata%i
    zs     =  simdata%z   !no that is not a typo

    ims = is(sim_ind:nYears, :)
    vms = vs(sim_ind:nYears, :)
    cms = cs(sim_ind:nYears, :)
    pims = zs(sim_ind:nYears, :)
    dms = ds(sim_ind:nYears, :)
    pifs = zs(sim_ind-1:nYears-1, :)
    cifs = cs(sim_ind-1:nYears-1, :)

    nobs = (nYears-sim_ind+1)*nFirms

    ! Look at some profitability values
    mu_z = sum(pims)/real(nobs)
    var_z = sum(pims**2.0_dp)/real(nobs)-mu_z**2.0_dp
    rho_z = sum( (pims - mu_z)*(pifs - mu_z) )/real(nobs-1)/var_z

    mu_i = sum(ims)/real(nobs)
    var_i = sum(ims**2.0_dp)/real(nobs)-mu_i**2.0_dp

    ! cash
    mu_c = sum(cms)/real(nobs)
    var_c = sum(cms**2.0_dp)/real(nobs)-mu_c**2.0_dp
    rho_c = sum( (cms - mu_c)*(cifs - mu_c) )/real(nobs-1)/var_c



    allmom(1)  = mu_c                             !allmoms(1)        Mean cash
    allmom(2)  = var_c                            !allmoms(2)        VAR cash
    allmom(3)  = mu_i                             !allmoms(3)        Mean investment
    allmom(4)  = var_i                            !allmoms(4)        VAR investment
    allmom(5)  = mu_z                             !allmoms(5)        Average profits
    allmom(6)  = var_z                            !allmoms(6)        VAR profits
    allmom(7)  = rho_z                            !allmoms(7)        autocorr profits


    momvector = allmom(pickout)

    deallocate(vs)
    deallocate(cs)
    deallocate(ds)
    deallocate(is)
    deallocate(pis)
    deallocate(zs)
    deallocate(vms)
    deallocate(cms)
    deallocate(dms)
    deallocate(ims)
    deallocate(pims)
    deallocate(pifs)
    deallocate(cifs)

end subroutine makemoments