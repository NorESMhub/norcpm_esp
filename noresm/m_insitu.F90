! File:          m_insitu.F90
!
! Created:       6 Feb 2008
!
! Last modified: 13 Feb 2008
!
! Author:        Pavel Sakov
!                NERSC
!
! Purpose:       The code to deal with insitu observations.
!
! Description:   This module contains the following subroutines:
!                  - insitu_setprofiles
!                      breaks the measurements into profiles and returns
!                      arrays of start and end indices for each profile
!                  - insitu_writeprofiles
!                      writes profiles to a netCDF file
!                  - insitu_prepareobs
!                      sorts out the measurements within profiles so they
!                      go in surface to bottom order and thins the measurements
!                      by keeping max 1 measurements per layer of the first
!                      ensemble member
!                It also contains the following data:
!                  nprof
!                    - the number of profiles
!                  pstart(nprof)
!                    - start indices for each profile in the array "obs" of
!                      type(measurement) stored in module m_obs
!                  pend(nprof)
!                    - end indices for each profile
!
! Modifications:
!                09/12/2024 Ping-Gin: For norcpm_esp, remove THIN block
!                18/03/2015 Yiguo Wang: set TEM_MIN = -2.5 and SAL_MIN = 1.0,
!                  which match with m_Generate_element_Si.F90
!                21/07/2015 Yiguo Wang: added a parameter (logical) THIN
!                30/7/2010 PS: added profile pivot points to profile output
!                  files (SAL.nc etc.)
!                29/7/2010 PS: some rather minor changes, including interface
!                  of insitu_writeforecast()
!                13/02/2008 PS: added insitu_writeprofiles()
!                26/02/2008 PS: put "nprof", "pstart" and "pend" as public data
!                  in this module
!                20/04/2008 PS: added insitu_QC() and insitu_writeforecast()
!                29/07/2010 PS: removed insitu_QC(). There is a generic obs QC
!                  procedure in m_obs.F90 now.

module m_insitu
  use mod_measurement
  use m_parse_blkdat
#if defined (QMPI)
   use qmpi
#else
   use qmpi_fake
#endif
  implicit none

  !
  ! public stuff
  !
  integer, allocatable, dimension(:), public :: pstart
  integer, allocatable, dimension(:), public :: pend
  integer, public :: nprof

  public insitu_setprofiles
  public insitu_prepareobs
  public insitu_writeprofiles

  !
  ! private stuff
  !

  real, parameter, private :: ONEMETER = 98060.0
  integer, parameter, private :: STRLEN = 512

  ! The portion of the layer thickness at which the variability in
  ! vertical data will be used for estimating the vertical representativeness
  ! error.
  !
  real, parameter, private :: VARCOEFF1 = 0.15
  
  ! A factor by which a calculated vertical representativeness error variance
  ! will be reduced if the data is in different layers
  !
  real, parameter, private :: VARCOEFF2 = 2.0

  ! Write information about this profile. Set to < 1 to switch off.
  !
  integer, parameter, private :: PDEBUGINFO = 0

  ! Integers used to tag the fields (to avoid parsing the string tags)
  !
  integer, parameter, private :: NONE = 0
  integer, parameter, private :: TEMPERATURE = 1
  integer, parameter, private :: SALINITY = 2

  real, parameter, private :: TEM_MIN = -2.5
  real, parameter, private :: TEM_MAX = 35.0
  real, parameter, private :: SAL_MIN = 1.0
  real, parameter, private :: SAL_MAX = 41.0

  ! Maximum allowed deviation between the observation and ensemble mean in
  ! terms of combined standard deviation.
  !
  real, parameter, private :: SAL_MAXRATIO = 10.0
  real, parameter, private :: TEM_MAXRATIO = 5.0

  ! If an observation is not considered an outlier, the observation error
  ! variance is modified so that the distance between the observation and the
  ! endemble mean is within DIST_MAX * sqrt(sigma_obs^2 + sigma_ens^2).
  ! Bigger values of DIST_MAX result in a more violent assimilation.
  !
  real, parameter, private :: DIST_MAX = 2.0

  ! Thin the profile 
  logical, parameter, private :: THIN = .FALSE.

contains

  ! Work out the number of profiles, each identified by "obs % i_orig_grid"
  ! and return start id of the first and the last obs in the profile in
  ! arrays "pstart" and "pend". "pstart" and "pend" are publicly available
  ! arrays stored in this module.
  !
  subroutine insitu_setprofiles(obstag, nobs, obs)
    character(*), intent(in) :: obstag
    integer, intent(in) :: nobs
    type(measurement), dimension(:), intent(inout) :: obs

    integer, allocatable, dimension(:) :: tmp, tmp1
    integer :: o, o1, o2, p, nobsp
    type(measurement), allocatable, dimension(:) :: tmpobs

    if (nobs == 0) then
       return
    end if

    if (allocated(pstart)) then
       deallocate(pstart)
       deallocate(pend)
    end if

    ! find the very first obs of the right kind
    !
    o1 = 1
    do while (trim(obs(o1) % id) /= trim(obstag) .and. o1 <= nobs)
       o1 = o1 + 1
    end do

    if (o1 > nobs) then
       return
    end if

    ! find the very last obs of the right kind
    !
    o2 = nobs
    do while (trim(obs(o2) % id) /= trim(obstag) .and. o2 >= 0)
       o2 = o2 - 1
    end do

    nprof = 1
    do o = o1+1, o2
       if (obs(o) % ipiv /= obs(o - 1) % ipiv .or.&
            obs(o) % jpiv /= obs(o - 1) % jpiv .or.&
            obs(o) % date /= obs(o - 1) % date) then
          nprof = nprof + 1
       end if
    end do

    allocate(pstart(nprof))
    allocate(pend(nprof))

    ! identify profiles
    !
    ! PS: This is a tricky cycle but it seems it is doing the job. Do not
    ! meddle with it.
    !
    pend = 0
    nprof = 1
    pstart(1) = o1
    do o = o1, o2
       ! find obs from the same profile
       !
       if (trim(obs(o) % id) == trim(obstag) .and.&
            ((obs(o) % i_orig_grid > 0 .and.&
            obs(o) % i_orig_grid == obs(pstart(nprof)) % i_orig_grid) .or.&
            (obs(o) % i_orig_grid <= 0 .and.&
            obs(o) % ipiv == obs(pstart(nprof)) % ipiv .and.&
            obs(o) % jpiv == obs(pstart(nprof)) % jpiv .and.&
            obs(o) % date == obs(pstart(nprof)) % date))) then
          pend(nprof) = o
          cycle
       end if

       if (trim(obs(o) % id) /= trim(obstag)) then
          print *, 'ERROR: insitu_setprofiles(): obs id does not match processed obs tag'
          stop
       end if

       ! if there were no obs of the right type in this profile yet,
       ! then pend(nprof) has not been set yet and therefore the condition
       ! below will yield "false"
       !
       if (pend(nprof) >= pstart(nprof)) then
          nprof = nprof + 1
       end if

       if (PDEBUGINFO > 0) then
          print *, '  DEBUG: new profile #', nprof, ', o =', o, ', id =', obs(o) % i_orig_grid
       end if
       pstart(nprof) = o
       pend(nprof) = o
    end do
    if (pend(nprof) < pstart(nprof)) then
       nprof = nprof - 1
    end if
    ! truncate "pstat" and "pend" to length "nprof"
    !
    allocate(tmp(nprof))
    tmp = pstart(1 : nprof)
    deallocate(pstart)
    allocate(pstart(nprof))
    pstart = tmp
    tmp = pend(1 : nprof)
    deallocate(pend)
    allocate(pend(nprof))
    pend = tmp
    deallocate(tmp)

    ! for glider data - sort observations in each profile by increasing depth
    !
    if (trim(obstag) == 'GSAL'.or. trim(obstag) == 'GTEM') then
       allocate(tmp(nobs))
       allocate(tmp1(nobs))
       allocate(tmpobs(nobs))
       do p = 1, nprof
          nobsp = pend(p) - pstart(p) + 1
          do o = 1, nobsp
             tmp(o) = o
          end do
          !
          ! (using procedure from pre_local_analysis())
          !
          call order(dble(nobsp), obs(pstart(p) : pend(p)) % depth,&
               dble(nobsp), tmp, tmp1)
          tmpobs(1 : nobsp) = obs(pstart(p) : pend(p))
          do o = 1, nobsp
             obs(pstart(p) + o - 1) = tmpobs(tmp1(o))
          end do
       end do
       deallocate(tmp, tmp1, tmpobs)
    end if
  end subroutine insitu_setprofiles


  ! 1. Sort out the obs within profiles so that they are stored in order of
  !    increasing depth.
  ! 2. Option: thin observations by keeping a single obs within a layer using the
  !    layers from the first ensemble member
  !
  subroutine insitu_prepareobs(obstag, nobs, obs)
    character(*), intent(in) :: obstag
    integer, intent(inout) :: nobs
    type(measurement), dimension(:), intent(inout) :: obs

    ! profiles
    !
    integer, allocatable, dimension(:) :: pnow
    integer :: nobs_max

    integer :: p, o
    type(measurement), allocatable, dimension(:) :: profile

    integer, allocatable, dimension(:) :: ipiv, jpiv
    real, allocatable, dimension(:) :: a1, a2, a3, a4
    real, allocatable, dimension(:) :: z1, z2

    integer :: nrev
    integer :: ndel
    integer :: oo
    real :: rdummy
    integer :: k, nk, ni, nj
    character(80) :: fname
    integer :: tlevel
    real, allocatable, dimension(:, :) :: dz2d
    real, dimension(2, 2) :: dz_cell
    real :: dz, zcentre
    integer :: best
    logical :: isrogue

    ! As we thin the measurements within each layer, it still may be a good
    ! idea to update the obs error variance if the variability within the layer
    ! is big enough. `dmin' and `dmax' give the min and max measured values
    ! within the layer.
    !
    real :: dmin, dmax
    real :: var1, var2

    integer :: nobsnew, nobs_thistype, nobs_othertype
    
    if (master) then
       print '(a, a, a)', '   insitu_prepareobs(', trim(obstag), '):'
       print '(a, i6)', '     total # of obs = ', nobs
    end if

    if (nobs == 0) then
       return
    end if

    call insitu_setprofiles(trim(obstag), nobs, obs)

    if (master) then
       print '(a, a, a, i6)', '     # of obs of type "', trim(obstag), '" = ',&
            sum(pend(1 : nprof) - pstart(1 : nprof)) + nprof
       print '(a, i4)', '     # of profiles = ', nprof
    end if

    ! find the maximal # of obs in a single profile
    !
    nobs_max = 0
    do p = 1, nprof
       nobs_max = max(nobs_max, pend(p) - pstart(p) + 1)
    end do

    if (master) then
       print '(a, i4)', '     max # of obs in a profile before thinning = ', nobs_max
    end if

    ! reverse the obs in profiles that go from bottom to surface
    !
    allocate(profile(nobs_max))
    nrev = 0
    do p = 1, nprof
       if (obs(pstart(p)) % depth > obs(pend(p)) % depth) then
          
          profile(1 : pend(p) - pstart(p) + 1) = obs(pstart(p) : pend(p))
          do o = 0, pend(p) - pstart(p)
             obs(pstart(p) + o) = profile(pend(p) - o)
          end do
          nrev = nrev + 1
       end if
    end do
    deallocate(profile)

    if (nrev > 0 .and. master) then
       print *, '  ', nrev, ' profile(s) reversed'
    end if

    ! check for rogue obs
    !
    ndel = 0
    do p = 1, nprof
       isrogue = .false. 
       do o = pstart(p) + 1, pend(p)

          ! shift the remaining obs in this profile one obs down
          !
          if (obs(o) % depth <= obs(o - 1) % depth) then
             isrogue = .true. 
             do oo = o + 1, pend(p)
                obs(oo - 1) = obs(oo)
             end do
             ndel = ndel + 1
             pend(p) = pend(p) - 1
          end if
       end do
       if (isrogue .and. master) then 
          print *, '  a rogue obs detected in profile # ', p 
       end if  
    end do

    if (ndel > 0 .and. master) then
       print *, '  ', ndel, 'rogue obs deleted'
    end if

  end subroutine insitu_prepareobs


  ! Write profiles to a NetCDF file
  !
  subroutine insitu_writeprofiles(fname, obstag, nobs, obs)
    use nfw_mod
    
    character(*), intent(in) :: fname
    character(*), intent(in) :: obstag
    integer, intent(inout) :: nobs
    type(measurement), dimension(:), intent(inout) :: obs

    ! profiles
    !
    integer :: p
    integer :: npoints, npoints_max

    ! I/O
    !
    integer :: ncid
    integer :: nprof_id(1), nk_id(1), dids(2)
    integer :: lat_id, lon_id, ipiv_id, jpiv_id, npoints_id, depth_id, v_id, variance_id
    character(STRLEN) :: varname

    real(8), allocatable, dimension(:, :) :: v

    if (.not. allocated(pstart)) then
       call insitu_setprofiles(trim(obstag), nobs, obs)
    end if

    call nfw_create(fname, nf_write, ncid)

    call nfw_def_dim(fname, ncid, 'nprof', nprof, nprof_id(1))
    call nfw_def_var(fname, ncid, 'lat', nf_double, 1, nprof_id, lat_id)
    call nfw_def_var(fname, ncid, 'lon', nf_double, 1, nprof_id, lon_id)
    call nfw_def_var(fname, ncid, 'ipiv', nf_int, 1, nprof_id, ipiv_id)
    call nfw_def_var(fname, ncid, 'jpiv', nf_int, 1, nprof_id, jpiv_id)
    call nfw_def_var(fname, ncid, 'npoints', nf_int, 1, nprof_id, npoints_id)
    npoints_max = maxval(pend - pstart) + 1
    call nfw_def_dim(fname, ncid, 'nk', npoints_max, nk_id(1))
    dids(1) = nk_id(1)
    dids(2) = nprof_id(1)
    call nfw_def_var(fname, ncid, 'depth', nf_double, 2, dids, depth_id)
    if (trim(obstag) == 'SAL' .or. trim(obstag) == 'GSAL') then
       varname = 'salt'
    else if (trim(obstag) == 'TEM' .or. trim(obstag) == 'GTEM') then
       varname = 'temp'
    else
       varname = trim(obstag)
    end if
    call nfw_def_var(fname, ncid, trim(varname), nf_double, 2, dids, v_id)
    call nfw_def_var(fname, ncid, 'variance', nf_double, 2, dids, variance_id)

    call nfw_enddef(fname, ncid)

    call nfw_put_var_double(fname, ncid, lat_id, obs(pstart) % lat)
    call nfw_put_var_double(fname, ncid, lon_id, obs(pstart) % lon)
    call nfw_put_var_int(fname, ncid, ipiv_id, obs(pstart) % ipiv)
    call nfw_put_var_int(fname, ncid, jpiv_id, obs(pstart) % jpiv)
    call nfw_put_var_int(fname, ncid, npoints_id, pend - pstart + 1)

    ! depth
    !
    allocate(v(npoints_max, nprof))
    v = -999.0
    do p = 1, nprof
       npoints = pend(p) - pstart(p) + 1
       v(1 : npoints, p) = obs(pstart(p) : pend(p)) % depth
    end do
    call nfw_put_var_double(fname, ncid, depth_id, v)
    
    ! data
    !
    v = -999.0
    do p = 1, nprof
       npoints = pend(p) - pstart(p) + 1
       v(1 : npoints, p) = obs(pstart(p) : pend(p)) % d
    end do
    call nfw_put_var_double(fname, ncid, v_id, v)
    
    ! data error variance
    !
    v = -999.0
    do p = 1, nprof
       npoints = pend(p) - pstart(p) + 1
       v(1 : npoints, p) = obs(pstart(p) : pend(p)) % var
    end do
    call nfw_put_var_double(fname, ncid, variance_id, v)

    call nfw_close(fname, ncid)

    deallocate(v)
    deallocate(pstart)
    deallocate(pend)
  end subroutine insitu_writeprofiles


  ! This subroutine appends the interpolated ensemble mean and the ensemble
  ! error variance to the assimilated profile data SAL.nc or TEM.nc. It also
  ! overwrites the observation error variance with latest values.
  !
  subroutine insitu_writeforecast(obstag, nobs, nrens, S, obs)
    use nfw_mod
    implicit none

    character(*), intent(in) :: obstag
    integer, intent(in) :: nobs
    integer, intent(in) :: nrens
    real, dimension(nobs, nrens), intent(in) :: S
    type(measurement), dimension(nobs), intent(inout) :: obs
    
    character(STRLEN) :: fname
    real, dimension(nobs) :: Smean, Svar
    integer :: i, p

    integer :: ncid
    integer :: dids(2)
    integer :: v_id, variance_id
    integer :: npoints_max, npoints
    real(8), allocatable, dimension(:, :) :: v

    ! need to set profiles for the given observation type
    !
    call insitu_setprofiles(obstag, nobs, obs)

    write(fname, '(a, ".nc")') trim(obstag)
    print *, 'Appending interpolated forecast for "', trim(obstag),&
         '" to "', trim(fname), '"'

    Smean = sum(S, DIM = 2) / nrens
    Svar = 0.0
    do i = 1, nobs
       Svar(i) = sum((S(i, :) - Smean(i)) ** 2)
    end do
    Svar = Svar / real(nrens - 1)

    call nfw_open(fname, nf_write, ncid)
    
    call nfw_inq_dimid(fname, ncid, 'nk', dids(1))
    call nfw_inq_dimid(fname, ncid, 'nprof', dids(2))

    call nfw_redef(fname, ncid)

    call nfw_def_var(fname, ncid, 'forecast', nf_double, 2, dids, v_id)
    call nfw_def_var(fname, ncid, 'forecast_variance', nf_double, 2, dids, variance_id)

    call nfw_enddef(fname, ncid)

    npoints_max = maxval(pend - pstart) + 1
    allocate(v(npoints_max, nprof))

    v = -999.0
    do p = 1, nprof
       npoints = pend(p) - pstart(p) + 1
       v(1 : npoints, p) = Smean(pstart(p) : pend(p))
    end do
    call nfw_put_var_double(fname, ncid, v_id, v)

    v = -999.0
    do p = 1, nprof
       npoints = pend(p) - pstart(p) + 1
       v(1 : npoints, p) = Svar(pstart(p) : pend(p))
    end do
    call nfw_put_var_double(fname, ncid, variance_id, v)

    ! update observation error variance
    !
    call nfw_redef(fname, ncid)
    call nfw_rename_var(fname, ncid, 'variance', 'variance_orig')
    call nfw_def_var(fname, ncid, 'variance', nf_double, 2, dids, variance_id)
    call nfw_enddef(fname, ncid)

    v = -999.0
    do p = 1, nprof
       npoints = pend(p) - pstart(p) + 1
       v(1 : npoints, p) = obs(pstart(p) : pend(p)) % var
    end do
    call nfw_put_var_double(fname, ncid, variance_id, v)

    call nfw_close(fname, ncid)

    deallocate(v)
  end subroutine insitu_writeforecast

end module m_insitu
