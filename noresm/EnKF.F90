!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a dummy module, for running ESP component in NorESM2.
!       [2023-10] fork by Ping-Gin Chiu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! File:          EnKF.F90
!
! Created:       ???
!
! Last modified: 20/04/2010
!
! Purpose:       Main program for EnKF analysis
!
! Description:   The workflow is as follows:
!                -- read model parameters
!                -- read obs
!                -- conduct necessary pre-processing of obs (superobing)
!                -- calculate ensemble observations
!                -- calculate X5
!                -- update the ensemble
!
! Modifications:
!                27/08/2024 PG: adjust for get BLOM data from memory and mpi
!                15/09/2014 YW: Coupling automatically layer thicknesses to preserve 
!                  the non-negativity of DP. The list of modified files is as follows:
!                  -- distribute.F90
!                  -- m_local_analysis.F90                              
!                28/10/2011 FC: The code is adapted to work with micom
!                20/9/2011 PS:
!                  Modified code to allow individual inflations for each of
!                  `NFIELD' fields updated in a batch - thanks to Ehouarn Simon
!                  for spotting this inconsistency
!                6/8/2010 PS:
!                  Small changes in calls to calc_X5() and update_fields() to
!                  reflect changes in interfaces.
!                6/7/2010 PS:
!                  Moved point output to a separate module m_point2nc.F90
!                25/5/2010 PS:
!                  Added inflation as a 4th command line argument
!                20/5/2010 PS:
!                  Set NFIELD = 4. This requires 4 GB per node in TOPAZ and
!                  "medium" memory model on Hexagon (a single allocation for a
!                   variable over 2GB)
!                20/4/2010 PS:
!                  Set NFIELD = 4. This will require 2 GB per node in TOPAZ.
!                  Thanks to Alok Gupta for hinting this possibility.
!                10/4/2010 PS:
!                  Moved variable `field' from real(8) to real(4);
!                  set NFIELD = 2.
!                Prior history:
!                  Not documented.

module m_EnKF

contains

subroutine EnKF()

#if defined(QMPI)
  use qmpi
#else
  use qmpi_fake
#endif
  use m_parameters
  use distribute
  use mod_measurement
  use m_get_micom_nrens
  use m_get_micom_grid
  use m_get_micom_dim
  use m_obs
  use m_local_analysis
  use m_prep_4_EnKF
  use m_set_random_seed2
  use m_get_micom_fld
  use m_put_micom_fld
  use mod_analysisfields
  use m_parse_blkdat
  use m_random
  use m_point2nc
  use netcdf
  use nfw_mod
  use norcpm_otf, only: gather_blom_1i1l,who_need_data_raise_hand &
                      , get_pe_start_end_iter, totalpe            &
                      , put_blom_1i1l, reset_get_pe_start_end_iter&
                      , set_uf_comm, free_uf_comm, save_nc_fld_update_fields &
                      , bcast_x5_otf_from_pe0
  use perf_mod, only: t_startf,t_stopf !!,t_initf,t_finalizef

  implicit none

  character(*), parameter :: ENKF_VERSION = "2.11"

  integer, external :: iargc

  ! NFIELD is the number of fields (x N) passed for the update during a call to
  ! update_fields(). In TOPAZ4 NFIELD = 2 if there is 1 GB of RAM per node, and
  ! NFIELD = 4 if there are 2 GB of RAM. Higher value of NFIELD reduces the
  ! number of times X5tmp.uf is read from disk, which is the main bottleneck
  ! for the analysis time right now.
  !
  !! In norcpm_otf, NFIELD can be less so that data can be spread to more processes
  !!     but cause error in update_fields()
  integer, parameter :: NFIELD = 53 

  character(512) :: options

  integer :: nrens
  real, allocatable, dimension(:,:) :: modlon, modlat, depths, readfld, readfld2
  real, allocatable, dimension(:,:) :: S ! ensemble observations HE
  real, allocatable, dimension(:)   :: d ! d - Hx

  integer k, m

  ! "New" variables used in the parallelization 
  integer, dimension(:,:), allocatable :: nlobs_array
  real(4), allocatable :: fld(:,:),dpfld(:,:),fld_ave(:),nb_ave(:)
  real(8) rtc, time0, time1, time2

  ! Additional fields
  character(len=3) :: cmem
  character(len=80) :: memfile
  integer :: fieldcounter

  character(100) :: text_string

  real :: rdummy
  integer :: idm, jdm, kdm
  integer :: i

  real :: mindx
  real :: meandx
  integer :: m1, m2, nfields
  real :: infls(NFIELD)
  logical :: isdp(NFIELD)

  !! norcpm_otf mod
  logical :: ismyiter             !! based on distrubte.F90
  logical :: is_scatter           !! scatter updated fld
  integer :: dstrank,srcrank      !! rank of ALLESPID
  integer :: nmyiter
  logical :: isalldone , debug=.false.
  integer :: infloopn
  integer :: t0,t1,clockrate     !! timing

#if defined(QMPI)
  call start_mpi()
#endif

  ! Read the characteristics of the assimilation to be carried out.

  !!not use here!if (iargc() /= 1) then
  !!not use here!   print *, 'Usage: EnKF <parameter file>'
  !!not use here!   print *, '       EnKF -h'
  !!not use here!   print *, 'Options:'
  !!not use here!   print *, '  -h -- describe parameter fie format'
  !!not use here!   call stop_mpi()
  !!not use here!else
  !!not use here!   call getarg(1, options)
  !!not use here!   if (trim(options) == "-h") then
  !!not use here!      call prm_describe()
  !!not use here!      call stop_mpi()
  !!not use here!   end if
  !!not use here!end if

  if (master) then
     print *
     print '(a, a)', ' EnKF version ', ENKF_VERSION
     print *
  end if

  call prm_read()
  call prm_print()

  ! get model dimensions
  !29/05/2015 Add reading of kdm
  call get_micom_dim(idm,jdm,kdm)
  if (master) then
     print *, 'read dimension idm,jdm,kdm :',idm,jdm,kdm
  end if
  allocate(modlon(idm,jdm))
  allocate(readfld(idm,jdm))
  allocate(readfld2(idm,jdm))
  allocate(modlat(idm,jdm))
  allocate(depths(idm,jdm))
  allocate(nlobs_array(idm, jdm))
  ! get model grid
  !
  !!  should it be parallel?
  call get_micom_grid(modlon, modlat, depths, mindx, meandx, idm, jdm)
  if (master) then
     print *,'MEAN grid size and min from scpx/scpy :',meandx,mindx
  end if

  ! set a variable random seed
  !
  !call set_random_seed3

  ! initialise point output
  !
  call p2nc_init

  !call t_startf('EnKF:initialisation')
  time0 = rtc()

  ! read measurements
  !
  call t_startf('EnKF:obs_readobs()')
  if (master) then
     print *, 'EnKF: reading observations'
  end if
  nobs = -1 !! make obs_readobs() read new linked obs 
  call system_clock(t0,clockrate)
  call obs_readobs
  call system_clock(t1)
  if(master)print'(a,f8.4)','EnKF timing, obs_readobs(): ',real(t1-t0)/real(clockrate)

  if (master) then
     print '(a, i6)', '   # of obs = ', nobs
     print '(a, a, a, e10.3, a, e10.3)', '   first obs = "', trim(obs(1) % id),&
          '", v = ', obs(1) % d, ', var = ', obs(1) % var
     print '(a, a, a, e10.3, a, e10.3)', '   last obs = "', trim(obs(nobs) % id),&
          '", v = ', obs(nobs) % d, ', var = ', obs(nobs) % var
  end if
  if (master) then
     print *
  end if
  call t_stopf('EnKF:obs_readobs()')

  ! read ensemble size and store in A
  !
  nrens = get_micom_nrens(idm, jdm)
  if (master) then
     print '(a, i5, a)', ' EnKF: ', nrens, ' ensemble members found'
  end if
  if (ENSSIZE > 0) then
     ENSSIZE = min(nrens, ENSSIZE)
  else
     ENSSIZE = nrens
  end if
  if (master) then
     print '(a, i4, a)', ' EnKF: ', ENSSIZE, ' ensemble members used'
  end if
  if (master) then
     print *
  end if

  ! PS - preprocess the obs using the information about the ensemble fields
  ! here (if necessary), before running prep_4_EnKF(). This is necessary e.g.
  ! for assimilating in-situ data because of the dynamic vertical geometry in
  ! HYCOM
  !
  call t_startf('EnKF:obs_prepareobs()')
  call system_clock(t0,clockrate)
  call obs_prepareobs
  call system_clock(t1)
  if(master)print'(a,f8.4)','EnKF timing, obs_prepareobs(): ',real(t1-t0)/real(clockrate)
  call t_stopf('EnKF:obs_prepareobs()')

  allocate(S(nobs, ENSSIZE), d(nobs))
  call t_startf('EnKF:prep_4_EnKF()')
  call system_clock(t0,clockrate)
  call prep_4_EnKF(ENSSIZE, d, S, depths, meandx / 1000.0, idm, jdm, kdm)
  call system_clock(t1)
  if(master)print'(a,f8.4)','EnKF timing, prep_4_EnKF(): ',real(t1-t0)/real(clockrate)
  if (master) then
     print *, 'EnKF: finished initialisation, time = ',  rtc() - time0
  end if
  call t_stopf('EnKF:prep_4_EnKF()')
  !call t_stopf('EnKF:initialisation')

  ! (no parallelization was required before this point)

  time1 = rtc()

  call t_startf('EnKF:calc_X5()')
  allocate(X5(ENSSIZE, ENSSIZE, idm))
  allocate(X5check(ENSSIZE, ENSSIZE, idm))
  call system_clock(t0,clockrate)
  call calc_X5(ENSSIZE, modlon, modlat, depths, mindx, meandx, d, S,&
       LOCRAD, RFACTOR2, nlobs_array, idm, jdm) !! output to tmpX5.uf
  deallocate(d, S, X5check)
  call system_clock(t1)
  if(master)print'(a,f8.4)','EnKF timing, calc_X5(): ',real(t1-t0)/real(clockrate)
  call system_clock(t0,clockrate)
  call bcast_x5_otf_from_pe0()
  call system_clock(t1)
  if(master)print'(a,f8.4)','EnKF timing, bcast_x5_otf_from_pe0(): ',real(t1-t0)/real(clockrate)
  if (master) then
     print *, 'EnKF: finished calculation of X5, time = ', rtc() - time0
  end if

  allocate(fld(idm * jdm, ENSSIZE * NFIELD))
  call t_stopf('EnKF:calc_X5()')

#if defined(QMPI)
  call barrier()
#endif

  ! get fieldnames and fieldlevels
  !
  call get_analysisfields()

  call distribute_iterations_field(numfields,fieldnames,fieldlevel)
#if defined(QMPI)
  call barrier() !KAL - just for "niceness" of output
#endif
  time2 = rtc()

  !! data assimilation with parallel process

  !! OTF method 2, 
  !!!   1. fill fld of all tasks, each needs NFIELD layers
  !!!   2. all tasks do update_fields()
  !!!   3. put fld data to OCN tasks
  !!!   4. fill fld of all tasks again, and so on
  
  infloopn = 0
  call reset_get_pe_start_end_iter()
  do  !! infinite loop
  !!  feed loop, layer_count
  !!    ask 1 task iteration start/continue and end/NFIELD, remember the last
  !!    ocn tasks, feed him, until end/NFIELD
  !!    loop though tasks, cycle
    infloopn = infloopn +1
    isalldone = .true.
    nmyiter = 0
    call t_startf('EnKF:gather model data')
    !print'(a,i0)','EnKF.F90:329, start read data, infloopn = ',infloopn
    do dstrank = 0, totalpe-1
        call get_pe_start_end_iter(dstrank,my_first_iteration,my_last_iteration&
                                    &,NFIELD,1,m1,m2,ismyiter) 
        !if(ismyiter)write(*,*)'EnKF.F90:318, mype,m1,m2 = ',qmpi_proc_num,m1,m2
        if (m1.eq.0 .and. m2.eq.0)cycle
        if (m1.gt.m2)cycle
        isalldone = .false.
        if(ismyiter) nmyiter = m2-m1+1
        do m = m1,m2  !! read fld
            if (ismyiter .and. trim(fieldnames(m)) /= 'dp' &
                         .and. fieldlevel(m)>=3 &
                         .and. .not. allocated(dpfld)) then
                allocate(dpfld(idm*jdm,ENSSIZE) &
                        ,fld_ave(idm*jdm),nb_ave(idm*jdm),source=0.0_sp)
            endif ! not dp and depth .ge.3, need be fill
            if(ismyiter)then !! for update_fields()
                infls(m - m1 + 1) = prm_getinfl(trim(fieldnames(m)))
                isdp(m - m1 + 1) = (trim(fieldnames(m)) == 'dp' ) 
            end if
            do k = 1, ENSSIZE  !! read fld for each instance
                ! reshaping and conversion to real(4)

                !print'(a,i0)','EnKF.F90:315, call gather_blom_1i1l(), infloopn = ',infloopn
                call gather_blom_1i1l(k,fieldnames(m),fieldlevel(m),dstrank,readfld)
                !print'(a,i0)','EnKF.F90:317, gather_blom_1i1l() done, infloopn = ',infloopn
                if(ismyiter) fld(:, ENSSIZE*(m-m1)+k) = reshape(readfld, (/idm*jdm/))
                if(ismyiter.and.debug) then
                    call save_nc_fld_update_fields(readfld,trim(fieldnames(m)),fieldlevel(m),'before_uf')
                end if
                if ( trim(fieldnames(m)) /= 'dp'.and.fieldlevel(m)>=3) then !! fill missing with ensmean
                      !Do not apply fix in the ML dp=1..2
                    call gather_blom_1i1l(k,'dp',fieldlevel(m),dstrank,readfld2)
                    if(ismyiter)then
                        ! reshaping and conversion to real(4)
                        dpfld(:, k) = reshape(readfld2, (/idm * jdm/))
                        !10 cm
                        where(dpfld(:, k)>9806.)
                            fld_ave=fld_ave+reshape(readfld, (/idm * jdm/))
                            nb_ave=nb_ave+1
                        endwhere
                    end if
                end if !! fill missing with ensmean
            end do !! do k = 1, ENSSIZE
            if (ismyiter .and. trim(fieldnames(m)) /= 'dp' .and. fieldlevel(m)>=3) then !! fill missing with ensmean
                do k = 1, ENSSIZE  !! fill with average
                    do i = 1, idm*jdm
                        !10 cm
                        if( dpfld(i, k)<9806. .and. nb_ave(i)>0 ) then
                           fld(i, ENSSIZE * (m - m1) + k)= fld_ave(i)/nb_ave(i);
                        endif
                     enddo !!  do i = 1, idm*jdm
                enddo !! do k = 1, ENSSIZE
                deallocate(dpfld,fld_ave,nb_ave)
            endif !! fill missing with ensmean
            !!call p2nc_storeforecast(idm, jdm, ENSSIZE, numfields, m, fld(:, ENSSIZE * (m - m1) + 1 : ENSSIZE * (m + 1 - m1)))
            if(ismyiter)then
                call p2nc_storeforecast(idm, jdm, ENSSIZE, nmyiter, m, fld(:,:nmyiter))
            end if
        end do !! do m = m1,m2
    end do ! dstrank
    !print'(a,i0)','EnKF.F90:390, data read, infloopn = ',infloopn
    call t_stopf('EnKF:gather model data')
    if (isalldone) exit !! exit from the infinite loop
  !!  every task do update_fields(), ok not every
    !print'(a,i0)','EnKF.F90:394, before update_fields(), infloopn = ',infloopn
    call t_startf('EnKF:update_fields()')
    if (nmyiter.gt.0)then  !! if there's data to update
        !! update_fields(): read tmpX5.uf and cal fld change
        call set_uf_comm(1) !! set ufcomm ufrank as group 1
        call update_fields(idm, jdm, ENSSIZE, nmyiter, nlobs_array, depths,&
                            &fld(1,1), infls, isdp)
    else
        call set_uf_comm(2) !! set ufcomm ufrank as group 2, not used here
    end if !! if there's data to update
    call free_uf_comm() !! free ufcomm
    call t_stopf('EnKF:update_fields()')
    !print'(a,i0)','EnKF.F90:406, after update_fields(), infloopn = ',infloopn
    !print'(a,i0)','EnKF.F90:354, data updated, putting, infloopn = ',infloopn
  !!  putdata loop
  !!    scatter data from 1st rank
  !!  if layer_count .eq. numfields, exit loop
    !print'(a,i0)','EnKF.F90:411, before put model data, infloopn = ',infloopn
    call t_startf('EnKF:put model data')
    do srcrank = 0, totalpe-1
        call get_pe_start_end_iter(srcrank,my_first_iteration,my_last_iteration&
                                    &,NFIELD,2,m1,m2,ismyiter)
        if (m1.eq.0 .and. m2.eq.0)then
            print*,'EnKF.F90:417, m1,m2 = 0, cycle'
            cycle
        end if
        do m = m1,m2
            do k = 1, ENSSIZE
                ! reshaping and conversion to real(8)
                readfld = reshape(fld(:, ENSSIZE * (m - m1) + k), (/idm, jdm/))
                if(ismyiter.and.debug) &
                    call save_nc_fld_update_fields(readfld,trim(fieldnames(m)),fieldlevel(m),'after_uf')
                !print'(a,x,a,x,i0,x,g0,x,3i3)','EnKF.F90:423, before put_blom_1i1l(), fldname, fldlev,maxval(readfld),m1,m2 :' &
                    !,trim(fieldnames(m)),fieldlevel(m),maxval(readfld),m1,m2,k
                call system_clock(t0,clockrate)
                call put_blom_1i1l(k,fieldnames(m),fieldlevel(m),srcrank,readfld)
                call system_clock(t1)
                !print'(a,3(i0,x),x,f8.4)','EnKF.F90:426, after put_blom_1i1l(), infloopn,m,k = ',infloopn,m,k,real(t1-t0)/real(clockrate)
                !call barrier()
            end do ! do k = 1, ENSSIZE
        end do ! m = m1,m2
    end do ! srcrank

    call t_stopf('EnKF:put model data')
    !print'(a,i0)','EnKF.F90:374, data putted, next loop, infloopn = ',infloopn
  end do !! infinite loop

  deallocate(X5)
  deallocate(fld)

  call p2nc_writeforecast

  ! Barrier only necessary for timings
#if defined(QMPI)
  call barrier()
#endif
  if (master) then
     print *, 'EnKF: time for initialization = ', time1 - time0
     print *, 'EnKF: time for X5 calculation = ', time2 - time1
     print *, 'EnKF: time for ensemble update = ', rtc() - time2
     print *, 'EnKF: total time = ', rtc() - time0
  end if
#if defined(QMPI)
  call barrier()
#endif
  !print *, 'EnKF: Finished'
  !call stop_mpi()

end subroutine EnKF

#if defined(_G95_)
! not tested! - PS
!
real function rtc()
  integer :: c

  call system_clock(count=c)
  rtc = dfloat(c)
end function rtc
#endif

end module m_EnKF
