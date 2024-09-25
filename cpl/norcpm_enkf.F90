!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module
!
! Structure:
!
! Revisions:
!   2023-11 Ping-Gin Chiu: File created
!   2024-09 PG: adjust for norcpm_otf and daily DA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module norcpm_enkf

    use m_EnKF, only: EnKF
    use m_prep_obs, only: p_prep_obs
    !use m_micom_ensemble_init, only: micom_ensemble_init
    use m_ensave, only: p_ensave
    use mpi
    use netcdf
    use shr_sys_mod, only : shr_sys_abort
    use ESMF, only: ESMF_Clock
    use seq_timemgr_mod, only: seq_timemgr_EClockGetData
    use perf_mod, only: t_startf,t_stopf !!,t_initf,t_finalizef
    use norcpm_otf, only: amiinocn,mxd1,mxd2,gather_blom_1i1l_allo,iulog &
                          ,reset_obs_record,sync_obs_records &
                          ,perturb_temp_with_random_seed,test_xca_put_get

    implicit none
    save
    private
    logical :: first = .true.
    integer :: mpicomm, mype, ninst
    integer :: mpisize
    integer :: rootpe = 0
    logical :: master = .false.
    logical :: run_micom_init = .false. !! disabled, too many var conflicts
    logical :: da_time_test = .true.
    !! perturb temp when start up, set in esp_init_mct()
    logical,public :: do_perturb_temp = .false.

    !! namelist
    character(len=255)             :: INPUTDATA,OCNGRIDFILE,MEAN_MOD_DIR,RES
    !!!! fortran cannot auto-allocate array via reading namelist
    !!!! need a better solution
    character(len=255),dimension(3) :: &  
        PRODUCERLIST, OBSLIST, REF_PERIODLIST, MONTHLY_ANOM, COMBINE_ASSIM, &
        FREQUENCYLIST
    integer :: ANOMALYASSIM, ANOM_CPL, OSAS, RFACTOR, ENSAVE, ENSSIZE
    character(len=4)   :: fforano
    character(len=255) :: ANALYSIS_DIRNAME,RESULT_DIRNAME
    integer :: DA_SINCE_DATE
    !! namelist modelio
    character(len=255) :: diri,diro,logfile

    public :: norcpm_assim_step_otf

    !! the order of variable in namelist is meanful, need be same as esp_in
    namelist /norcpm/INPUTDATA,OCNGRIDFILE &
        ,PRODUCERLIST, OBSLIST,FREQUENCYLIST, REF_PERIODLIST, MONTHLY_ANOM &
        , COMBINE_ASSIM ,ANOMALYASSIM, ANOM_CPL, OSAS, RFACTOR,fforano & 
        , ANALYSIS_DIRNAME,RESULT_DIRNAME,RES,MEAN_MOD_DIR,ENSAVE,ENSSIZE &
        , DA_SINCE_DATE
    !! modelio is not imply yet
    namelist /modelio/diri,diro,logfile

contains
    subroutine init_save_variables()
        integer :: ierr
        integer :: seed
        real    :: perturb_amp

        call MPI_COMM_RANK(mpicomm,mype,ierr)
        call MPI_COMM_SIZE(mpicomm,mpisize,ierr)
        call if_true_then_msg_stop(ierr.ne.MPI_SUCCESS, &
                    'norcpm ERROR: MPI_COMM_RANK() error')
        if(mype .eq. rootpe)master = .true.

        call read_namelist_norcpm('esp_in')

        if(master) call EXECUTE_COMMAND_LINE('rm -rf   '//trim(ANALYSIS_DIRNAME))
        if(master) call EXECUTE_COMMAND_LINE('mkdir -p '//trim(ANALYSIS_DIRNAME))

        !! perturb temperature if decieded in esp_comp_mct.F90
        seed = mype !! same seed produce same perturbation
        perturb_amp = 1.e-8 !! max amplitude
        if(do_perturb_temp) then
            if(master)write(iulog,*)'  NorCPM: perturb temp: ' &
                ,perturb_amp
            call perturb_temp_with_random_seed(seed,perturb_amp)
            do_perturb_temp = .false.
        end if

        first = .false.
    end subroutine init_save_variables


    subroutine norcpm_assim_step_otf(mpicomm,EClock)
        !! this subroutine called by esp_run_mct()
        integer, intent(in) :: mpicomm !! all ESP MPI comm
        type(ESMF_Clock),intent(in) :: EClock
        logical :: should_return = .false.
        integer :: ierr, EnKF_CNT = 0,fid
        character(len=2) :: ec_str !! str for EnKF_CNT
        integer :: year,month,day,sec, ymd
        integer :: i
        logical :: read_dp, doDA
        character(len=16) :: datestr

        if (first) call init_save_variables()

        !! test test_xca_put_get(dstinst)
        !do i = 1,4
        !    if(mype.eq.0)print*,'norcpm_assim_step_otf(): test_xca_put_get()',i
        !    call test_xca_put_get(i)
        !end do
        !! get model date time
        call seq_timemgr_EClockGetData(EClock, curr_yr=year, curr_mon=month,&
            curr_day=day,curr_tod=sec,curr_ymd=ymd)

        if(ymd .lt. DA_SINCE_DATE)return !! before date to DA
        if(.not.is_any_da_time(FREQUENCYLIST,year,month,day,sec))return

        !! if matched any DA frequency
        write(datestr,'(i4.4,2(a,i2.2),a,i5.5)')year,'-',month,'-',day,'-',sec
        if(master) write(iulog,'(3a)',advance='no')'  NorCPM DA ',datestr,':'

        call MPI_BARRIER(mpicomm,ierr)

        !! prep_obs*3, distribute tasks
        EnKF_CNT = 0

        !! read obs data
        call reset_obs_record()
        read_dp = .true.
        doDA = .false.
        do i = 1, size(OBSLIST) 
            !! check run this da or not 
            if(is_da_time(FREQUENCYLIST(i),year,month,day,sec))then
                doDA = .true.
                if(master)write(iulog,'(2a)',advance='no')' ',trim(OBSLIST(i))

                if(master)print'(a,i0)','norcpm_enkf.F90:98, prep_obs, i=',i
                !! link data files
                !!!! read dp for read_EN4_profile()
                if(PRODUCERLIST(i)(1:3).eq.'EN4' .and. read_dp)then
                    call gather_blom_1i1l_allo(1,'dp',1,rootpe,mxd1)
                    call gather_blom_1i1l_allo(1,'dp',2,rootpe,mxd2)
                    read_dp = .false.
                end if
                call MPI_BARRIER(mpicomm,ierr)
                !! do links
                if(trim(FREQUENCYLIST(i)).eq.'DAY')then
                    call prep_prep_obs(trim(OBSLIST(i)),trim(PRODUCERLIST(i)) &
                        ,trim(REF_PERIODLIST(i)),trim(COMBINE_ASSIM(i))       &
                        ,year,month,day)
                else
                    call prep_prep_obs(trim(OBSLIST(i)),trim(PRODUCERLIST(i)) &
                        ,trim(REF_PERIODLIST(i)),trim(COMBINE_ASSIM(i))       &
                        ,year,month)
                end if
                call MPI_BARRIER(mpicomm,ierr) !! wait for small files prepared
                if(master)then
                    call exec_in_analysis_dir('prep_obs')
                    !call exec_in_analysis_dir(  &
                    !    'mv observations.uf observations.uf_' & 
                    !    //trim(OBSLIST(i))//'.'//trim(PRODUCERLIST(i)))
                end if
                call MPI_BARRIER(mpicomm,ierr)
            end if !(is_da_time(FREQUENCYLIST(i),year,month,day,sec))
        end do ! i = 1, size(OBSLIST) 

        !!not use, remove?!! if (trim(COMBINE_ASSIM(i)).ne.'1')cycle
        if (.not.doDA)return
        call sync_obs_records() !! bcast obs

        EnKF_CNT = EnKF_CNT+1
        if(master)then
            !write(iulog,'(a,i2)')'norcpm: EnKF_CNT= ',EnKF_CNT
            !flush(iulog)
            write(ec_str,'(i2)')EnKF_CNT
            ec_str = adjustl(ec_str)
            call write_analysisfields_file(trim(ANALYSIS_DIRNAME)//'/analysisfields.in')
            call write_enkf_prm_file(ENSSIZE,RFACTOR,trim(ANALYSIS_DIRNAME)//'/enkf.prm')
            !call exec_in_analysis_dir( &
            !    & 'cat observations.uf_* > observations.uf ; &
            !    &  rm -f observations.uf_*')
        end if ! master
        if(.false..and. ENSAVE.eq.1)then !! no ensave for now
            call exec_in_analysis_dir('ensave','forecast')
            if(master .and. ENSAVE.eq.1)then
                call exec_in_analysis_dir('mv forecast_avg.nc forecast_avg_'//trim(ec_str)//'.nc')
                ! no aicen in BLOM, skip ensave_ice()
            end if
        end if
        call MPI_BARRIER(mpicomm,ierr)
        if(master)print*,'norcpm_enkf.F90: 130, start EnKF'
        call exec_in_analysis_dir('EnKF')
        if(master)print*,'norcpm_enkf.F90: 132, EnKF done'
        call MPI_BARRIER(mpicomm,ierr)
        if(master)then
            call exec_in_analysis_dir( '&
                & mv enkf_diag.nc enkf_diag_'//trim(ec_str)//'.nc; &
                & mv tmpX5.uf tmpX5_'//trim(ec_str)//'.uf ')
        end if
        if(.false..and. ENSAVE.eq.1)then !! no ensave in otf for now
            call exec_in_analysis_dir('ensave','forecast')
            if(master) &
                call exec_in_analysis_dir('mv forecast_avg.nc analysis_avg_'//trim(ec_str)//'.nc')
            ! no aicen in BLOM, skip ensave_ice()
        end if
        !if(master)write(iulog,'(a)')'norcpm: Finished with EnKF; call number :'//trim(ec_str)

        if(master) write(iulog,'(a)')'.'

        if(run_micom_init)then 
            call MPI_BARRIER(mpicomm,ierr)
            call exec_in_analysis_dir('micom_ensemble_init')
            call MPI_BARRIER(mpicomm,ierr)
        end if

        ! archive assimilation files, need revise
        if(master) call archive_assimilation(year,month)
    end subroutine norcpm_assim_step_otf

    subroutine archive_assimilation(year,month)
        integer, intent(in) :: year,month
        character(len=4) :: yyyy
        character(len=2) :: mm
        write(yyyy,'(i4.4)')year
        write(mm,  '(i2.2)')month
        call exec_in_analysis_dir('mkdir -p '//trim(RESULT_DIRNAME)//'/'//yyyy//'_'//mm)
        call exec_in_analysis_dir('mv enkf_diag_*.nc observations-*.nc tmpX5_*.uf '//trim(RESULT_DIRNAME)//'/'//yyyy//'_'//mm)
        if(ENSAVE.eq.1)then
            !call exec_in_analysis_dir('mv analysis_*avg_*.nc forecast_*avg_*.nc '//trim(RESULT_DIRNAME)//'/'//yyyy//'_'//mm)
            !call exec_in_analysis_dir('mv SAL.nc TEM.nc '//trim(RESULT_DIRNAME)//'/'//yyyy//'_'//mm)
        end if
        !! no aicen and vicen in BLOM restart, need revise
    end subroutine archive_assimilation

    subroutine exec_in_analysis_dir(cmd,arg1,arg2,noerr_)
        character(len=*),intent(in) :: cmd
        character(len=*),intent(in),optional :: arg1,arg2
        logical, intent(in), optional :: noerr_
        character(:), allocatable :: precmd
        character(len=255) :: rundir
        integer :: ierr
        logical :: noerr !! ignore error
        noerr = .true.
        if(present(noerr_))noerr = noerr_

        if (master)write(*,'(2a)')'norcpm exec_in_analysis_dir(): ',cmd
        !! empty line in ifort!! call getcwd(rundir,ierr)
        call get_environment_variable(name='PWD',value=rundir,status=ierr)
        if (ierr.ne.0)write(iulog,'(a,i5)')'norcpm ERROR, chdir() error: ',ierr
        if(trim(cmd).eq.'EnKF')then
            call t_startf('EnKF')
            call chdir(trim(ANALYSIS_DIRNAME))
            call EnKF()
            call chdir(rundir)
            call t_stopf('EnKF')
            return
        end if
        if(trim(cmd).eq.'ensave')then
            call t_startf('ensave')
            call chdir(trim(ANALYSIS_DIRNAME))
            call p_ensave(trim(arg1),ENSSIZE)  !! need to change
            call chdir(rundir)
            call t_stopf('ensave')
            return
        end if
        if(trim(cmd).eq.'prep_obs')then
            call t_startf('prep_obs')
            !!write(iulog,'(2a)')'norcpm: chdir to: ',trim(ANALYSIS_DIRNAME)
            call chdir(trim(ANALYSIS_DIRNAME))
            call p_prep_obs() !! cannot be parallel, depend on i/o file names
                            !! need fill arguments
            call chdir(rundir)
            call t_stopf('prep_obs')
            return
        end if

        if(trim(cmd).eq.'micom_ensemble_init')then !! not done yet, too many conflict
            !if(amiinocn)then !! only run on OCN pes
            !    call t_startf('micom_ensemble_init')
            !    call chdir(trim(ANALYSIS_DIRNAME))
            !    call micom_ensemble_init(ENSSIZE)
            !    call chdir(rundir)
            !    call t_stopf('micom_ensemble_init')
            !end if
            return
        end if

        precmd = 'cd '//trim(ANALYSIS_DIRNAME)//'/ ;'
        call EXECUTE_COMMAND_LINE(precmd//' '//cmd, exitstat=ierr)
        if(ierr.ne.0.and. .not.noerr) &
            call shr_sys_abort('norcpm ERROR: cmd: '//precmd//' '//cmd)
    end subroutine exec_in_analysis_dir

    subroutine prep_prep_obs(OBSTYPE,PRODUCER,REF_PERIOD,COMB_ASSIM,year,month,day,sec)
        character(len=*),intent(in) :: OBSTYPE,PRODUCER,REF_PERIOD,COMB_ASSIM
        integer, intent(in)         :: year,month
        integer, intent(in),optional:: day,sec
        character(len=255) :: fn,fn2,dst,msg,cmd
        character(len=2) :: mm,dd
        character(len=4) :: yyyy
        character(len=5) :: sssss !! sec of day
        character(len=:),allocatable :: datestr
        write(mm,'(i2.2)')month
        write(yyyy,'(i4.4)')year
        if(present(day))write(dd,'(i2.2)')day
        if(present(sec))write(sssss,'(i5.5)')sec

        if(mype.eq.0) then ! obs data
            fn = 'INPUTDATA/obs/OBSTYPE/PRODUCER/YYYY_MM_DD.nc'
            dst= 'ANALYSIS_DIRNAME/YYYY_MM_DD.nc'
            fn = Replace_Text(fn,'INPUTDATA',trim(INPUTDATA))
            fn = Replace_Text(fn,'OBSTYPE'  ,trim(OBSTYPE))
            fn = Replace_Text(fn,'PRODUCER' ,trim(PRODUCER))
            fn = Replace_Text(fn,'YYYY'     ,yyyy)
            fn = Replace_Text(fn,'MM'       ,mm)
            if(present(day)) fn = Replace_Text(fn,'DD'       ,dd)
            fn = Replace_Text(fn,'_DD.nc'   ,'.nc') !! if no DD present

            dst= Replace_Text(dst,'ANALYSIS_DIRNAME',ANALYSIS_DIRNAME)
            dst= Replace_Text(dst,'YYYY'     ,yyyy)
            dst= Replace_Text(dst,'MM'       ,mm)
            if(present(day)) dst = Replace_Text(dst,'DD'       ,dd)
            dst= Replace_Text(dst,'_DD.nc'   ,'.nc') !! if no DD present
            msg= 'norcpm ERROR: obs data missing: '//trim(fn)
            call link_if_exist_or_stop(fn=trim(fn),dst=trim(dst),stopmsg=trim(msg))
        end if

        if(mype.eq.1) then ! model mean, use monthly if diaily is not present
            fn = 'MEAN_MOD_DIR/Free-averageMMdDD-REF_PERIOD.nc'
            fn = Replace_Text(fn,'MEAN_MOD_DIR',trim(MEAN_MOD_DIR))
            fn = Replace_Text(fn,'MM'     ,mm)
            fn = Replace_Text(fn,'REF_PERIOD',REF_PERIOD)

            fn2= Replace_Text(fn,'dDD-'     ,'-') !! if no DD present
            if(present(day))fn = Replace_Text(fn,'DD'     ,dd)
            fn = Replace_Text(fn,'dDD-'     ,'-') !! if no DD present
            msg= 'norcpm ERROR: model mean data missing: '//trim(fn)
            dst= trim(ANALYSIS_DIRNAME)//'/mean_mod.nc'
            call link_if_exist_or_stop(fn=trim(fn),fn2=trim(fn2),dst=trim(dst),stopmsg=trim(msg))
        end if

        if(mype.eq.2) then ! obs unc
            fn = 'INPUTDATA/enkf/RES/PRODUCER/RES_OBSTYPE_obs_unc_FFORANO.nc'
            fn = Replace_Text(fn,'INPUTDATA',trim(INPUTDATA))
            fn = Replace_Text(fn,'RES'      ,trim(RES))
            fn = Replace_Text(fn,'PRODUCER' ,trim(PRODUCER))
            fn = Replace_Text(fn,'OBSTYPE'  ,trim(OBSTYPE))
            fn = Replace_Text(fn,'FFORANO'  ,trim(fforano))

            !fn = trim(INPUTDATA)//'/enkf/'//trim(RES)//'/'//PRODUCER//'/'//&
            !    &trim(RES)//'_'//OBSTYPE//'_obs_unc_'//fforano//'.nc'

            dst= trim(ANALYSIS_DIRNAME)//'/obs_unc_'//OBSTYPE//'.nc'
            msg= 'norcpm: obs unc missing, but maybe fine: '//trim(fn)
            call link_if_exist_or_stop(fn=trim(fn),dst=trim(dst),stopmsg=trim(msg),stoprun=.false.)
        end if

        if(mype.eq.3) then ! obs mean
            datestr = mm
            if(present(day)) datestr = mm//'d'//dd
            fn ='INPUTDATA/obs/OBSTYPE/PRODUCER/OBSTYPE_avg_MMdDD-REF_PERIOD.nc'
            fn = Replace_Text(fn,'INPUTDATA',trim(INPUTDATA))
            fn = Replace_Text(fn,'OBSTYPE'  ,trim(OBSTYPE))
            fn = Replace_Text(fn,'PRODUCER' ,trim(PRODUCER))
            fn = Replace_Text(fn,'YYYY'     ,yyyy)
            fn = Replace_Text(fn,'MM'       ,mm)
            if(present(day))fn = Replace_Text(fn,'DD'     ,dd)
            fn = Replace_Text(fn,'dDD-'     ,'-') !! if no DD present
            fn = Replace_Text(fn,'REF_PERIOD',REF_PERIOD)
            !fn = trim(INPUTDATA)//'/obs/'//OBSTYPE//'/'//PRODUCER//'/'//OBSTYPE//'_avg_'//datestr//'-'//REF_PERIOD//'.nc'
            dst= trim(ANALYSIS_DIRNAME)//'/mean_obs.nc'
            msg= 'norcpm ERROR: obs mean missing: '//trim(fn)
            call link_if_exist_or_stop(fn=trim(fn),dst=trim(dst),stopmsg=trim(msg))
        end if

        if(mype.eq.4) then ! ocean grid file
            fn = trim(OCNGRIDFILE)
            dst = trim(ANALYSIS_DIRNAME)//'/grid.nc'
            msg= 'norcpm ERROR: grid file missing: '//trim(fn)
            call link_if_exist_or_stop(fn=trim(fn),dst=trim(dst),stopmsg=trim(msg))
        end if

        if(mype.eq.5) then ! obs data contain
            datestr = mm
            if(present(day)) datestr = mm//'_'//dd
            fn = 'INPUTDATA/enkf/infile.data.OBSTYPE.PRODUCER'
            fn = Replace_Text(fn,'INPUTDATA',trim(INPUTDATA))
            fn = Replace_Text(fn,'OBSTYPE'  ,trim(OBSTYPE))
            fn = Replace_Text(fn,'PRODUCER' ,trim(PRODUCER))
            !fn = trim(INPUTDATA)//'/enkf/infile.data.'//trim(OBSTYPE)//'.'//trim(PRODUCER)
            dst = trim(ANALYSIS_DIRNAME)//'/infile.data'
            cmd = 'sed -e "s/yyyy/'//yyyy//'/" -e "s/mm/'//datestr//'/" "'//trim(fn)//'" > "'//trim(dst)//'"'
            write(*,'(a)')'norcpm: run cmd='//trim(cmd)
            call EXECUTE_COMMAND_LINE(trim(cmd))
        end if
    end subroutine prep_prep_obs

    subroutine link_if_exist_or_stop(fn,dst,stopmsg,fn2,stoprun)
        character(len=*), intent(in) :: fn,dst
        character(len=*), intent(in),optional :: stopmsg,fn2
        character(len=255) :: pwd 
        logical, intent(in),optional :: stoprun
        logical :: stoprun_
        integer :: ierr = 0

        stoprun_ = .true.
        if(present(stoprun))stoprun_ = stoprun
        write(*,'(a)')'norcpm: soft link: '//fn//' -> '//dst
        if(is_file_exist(fn))then
            write(*,'(a)')"ln -sfr "//fn//' '//dst
            call EXECUTE_COMMAND_LINE("ln -sfr "//fn//' '//dst,exitstat=ierr)
            if (ierr.eq.0)return
        end if
        if(present(fn2))then
            if(is_file_exist(fn2))then
                write(*,'(a)')"ln -sfr "//fn2//' '//dst
                call EXECUTE_COMMAND_LINE("ln -sfr "//fn2//' '//dst,exitstat=ierr)
                if (ierr.eq.0)return
            end if
        end if
        if (ierr.ne.0)then
            call get_environment_variable(name='PWD',value=pwd,status=ierr)
            write(*,'(2a)')'norcpm ERROR: link process error,pwd=',trim(pwd)
            write(*,'(2a)')'norcpm ERROR: link process error,src=',trim(fn)
            write(*,'(2a)')'norcpm ERROR: link process error,dst=',trim(dst)
        end if
        if(stoprun_.and.present(stopmsg)) call shr_sys_abort(stopmsg)
        if(stoprun_) call shr_sys_abort('norcpm ERROR: link_if_exist_or_stop()')
    end subroutine link_if_exist_or_stop

    logical function is_file_exist(fn)
        character(len=*),intent(in) :: fn
        inquire(file=fn,exist=is_file_exist)
    end function is_file_exist

    logical function is_da_time(dafreq,year,month,day,sec)
        !! check is the model time to do DA
        !! cannot do DA when sec = 0, it cause esp run after restart output
        character(len=*), intent(in) :: dafreq
        integer,          intent(in) :: year,month,day,sec
        !local
        integer :: midday = 43200 !! set DA at middle of day
        is_da_time = .false.
        
        SELECT CASE(trim(dafreq))
            CASE('MONTH') !! 
                if(day.eq.15.and.sec.eq.midday)is_da_time = .true.
            CASE('3DAY') !! 
                if(mod(day,3).eq.0.and.sec.eq.midday)is_da_time = .true.
            CASE('DAY')
                if(sec.eq.midday)is_da_time = .true.
            CASE DEFAULT
                print*,'ERROR, is_da_time(): dafreq invaild: ',trim(dafreq)
                stop
        END SELECT
    end function is_da_time
    logical function is_any_da_time(dafreqs,year,month,day,sec)
        character(len=*), intent(in) :: dafreqs(:)
        integer,          intent(in) :: year,month,day,sec
        integer :: i
        do i = 1,size(dafreqs)
            if(is_da_time(dafreqs(i),year,month,day,sec)) then
                is_any_da_time = .true.
                return
            end if
        end do
        is_any_da_time = .false.
    end function is_any_da_time

    subroutine delink_and_link(src,dst)
        character(len=*), intent(in) :: src,dst
        character(len=255) :: cmd
        integer :: exitstat
        cmd = "test $(stat -c '%h' '"//src//"') -eq 1 && test ! -h '"//src//"'"
        write(*,'(a)')trim(cmd)
        call EXECUTE_COMMAND_LINE(cmd,exitstat=exitstat)
        if(exitstat.ne.0)then !! if src reference to other file
            call EXECUTE_COMMAND_LINE('cp "'//src//'" "'//src//'_dst"')
            call EXECUTE_COMMAND_LINE('mv --backup=t "'//src//'" "'//src//'_orig"')
            call EXECUTE_COMMAND_LINE('mv "'//src//'_dst" "'//src//'"')
        end if

        call EXECUTE_COMMAND_LINE("ln -sfr "//src//' '//dst)
    end subroutine delink_and_link

    subroutine read_namelist_norcpm(fn)
        character(len=*) :: fn
        integer :: fid, ios
        if (.true.)then !! the namelist generate is not ready yet, check buildnml
            open(file=fn,action='read',newunit=fid,iostat=ios)
            if(ios.ne.0) write(*,*)"norcpm, Error: "//trim(fn)//' namelist open error'
            read(nml=norcpm,iostat=ios,unit=fid)
            if(ios.ne.0) write(*,*)"norcpm, Error: "//trim(fn)//' namelist read error'
            close(fid)
        end if

        write(iulog,*)'norcpm: ENSSIZE = ',ENSSIZE
        ninst = ENSSIZE
        !open(file='esp_modelio.nml_0001',action='read',newunit=fid,iostat=ios)
        !read(nml=modelio,iostat=ios,unit=fid)
        !if(ios.ne.0)then
        !    write(iulog,*)"norcpm, Error: esp_modelio.nml_0001 namelist read error"
        !end if
        !close(fid)
        !open(file=logfile,action='write',form='formatted',newunit=iulog)
    end subroutine read_namelist_norcpm

    subroutine write_enkf_prm_file(ENSSIZE,RFACTOR,ofn)
        character(len=*), intent(in) :: ofn
        integer, intent(in) :: ENSSIZE,RFACTOR
        integer :: fid
        !! return string for enkf.prm
        !! need revise
        write(*,*)'norcpm: write_enkf_prm_file(): ENSSIZE = ',ENSSIZE
        open(file=ofn,newunit=fid)
        write(fid,'(a)')"&method"
        write(fid,'(a)')"  methodtag = 'DEnKF'"
        write(fid,'(a)')"/"
        write(fid,'(a)')"&ensemble"
        write(fid,'(a,i4)')"    enssize = ",ENSSIZE
        write(fid,'(a)')"/"
        write(fid,'(a)')"&localisation"
        write(fid,'(a)')"    locfuntag = 'Gaspari-Cohn'"
        write(fid,'(a)')"    locrad = 1500.0"
        write(fid,'(a)')"/"
        write(fid,'(a)')"&moderation"
        write(fid,'(a)')"    infl = 1.00"
        write(fid,'(a,i4)')"    rfactor1 = ",RFACTOR
        write(fid,'(a)')"    rfactor2 = 4.0"
        write(fid,'(a)')"    kfactor = 2.0"
        write(fid,'(a)')"/"
        write(fid,'(a)')"&files"
        write(fid,'(a)')"/"
        write(fid,'(a)')"&prmest /"
        write(fid,'(a)')"/"
        close(fid)
    end subroutine write_enkf_prm_file

    subroutine write_analysisfields_file(ofn)
        character(len=*),intent(in)  :: ofn
        integer :: fid
        open(file=ofn,newunit=fid)
        !write(fid,'(a)')'u         1 53'
        !write(fid,'(a)')'v         1 53'
        write(fid,'(a)')'dp        1 53'
        write(fid,'(a)')'temp      1 53'
        write(fid,'(a)')'saln      1 53'
        !!write(fid,'(a)')'uflx      1 53'
        !!write(fid,'(a)')'vflx      1 53'
        !!write(fid,'(a)')'utflx     1 53'
        !!write(fid,'(a)')'vtflx     1 53'
        !!write(fid,'(a)')'usflx     1 53'
        !!write(fid,'(a)')'vsflx     1 53'
        !write(fid,'(a)')'pb        1 1'
        !!write(fid,'(a)')'ub        1 1'
        !!write(fid,'(a)')'vb        1 1'
        !!write(fid,'(a)')'ubflx     1 1'
        !!write(fid,'(a)')'vbflx     1 1'
        !!write(fid,'(a)')'ubflxs    1 1'
        !!write(fid,'(a)')'vbflxs    1 1'
        !!write(fid,'(a)')'ubcors_p  0 0'
        !!write(fid,'(a)')'vbcors_p  0 0'
        !!write(fid,'(a)')'phi       0 0'
        !!write(fid,'(a)')'sealv     0 0'
        !!write(fid,'(a)')'ustar     0 0'
        !!write(fid,'(a)')'buoyfl    0 0 '
        close(fid)
    end subroutine write_analysisfields_file

    FUNCTION Replace_Text (s,text,rep)  RESULT(outs)
        !! from https://fortranwiki.org/fortran/show/String_Functions
        CHARACTER(*)        :: s,text,rep
        CHARACTER(LEN(s)+100) :: outs   ! provide outs with extra 100 char len
        !!CHARACTER,allocatable :: outs
        INTEGER             :: i, nt, nr

        outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
        DO
            i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
            outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
        END DO
    END FUNCTION Replace_Text


    logical function is_restarts_EnKF_already(ncfn)
        character(len=*), intent(in) :: ncfn
        is_restarts_EnKF_already = .false.
        write(*,*)'norcpm is restart file EnKFed? ',ncfn
        if (is_att_present(ncfn,NF90_GLOBAL,'enkf')) &
            is_restarts_EnKF_already = .true.
        return
    end function is_restarts_EnKF_already

    logical function is_att_present(ncfn,varid,att)
        character(len=*),intent(in) :: ncfn,att
        integer :: varid
        integer :: ierr, ncid
        !if(master)then
            call check_nc_err(nf90_open(ncfn,NF90_NOWRITE,ncid))
            ierr = NF90_INQUIRE_ATTRIBUTE(ncid,varid,trim(att))
            call check_nc_err(NF90_CLOSE(ncid))
            is_att_present = .false.
            if(ierr.eq.NF90_NOERR)is_att_present = .true.
        !end if 
        !call MPI_BCAST(is_att_present,1,MPI_LOGICAL,rootpe,mpicomm,ierr)
        return
    end function is_att_present

    subroutine check_nc_err(ierr,ierrstop)
        integer, intent(in) :: ierr
        logical, intent(in), optional :: ierrstop
        logical :: errstop
        errstop = .true.
        if(present(ierrstop))errstop = ierrstop
        if (ierr .ne. NF90_NOERR)then
            write(*,*)'norcpm ERROR: '//trim(NF90_STRERROR(ierr))
            if (errstop)call shr_sys_abort('norcpm ERROR: '//trim(NF90_STRERROR(ierr)))
        end if
        return
    end subroutine check_nc_err

    function current_restart_fn(comp,inst) result(fn)
        character(:), allocatable :: fn
        character(len=*), intent(in) :: comp
        character(len=255) :: tmp
        integer, intent(in) :: inst
        integer :: fid, ierr
        character(len=30) :: rpointer_fn

        write(rpointer_fn,'(a,i4.4)')'rpointer.'//trim(comp)//'_',inst
        open(file=rpointer_fn,action='read',newunit=fid,iostat=ierr)
        call if_true_then_msg_stop(ierr.ne.0,'norcpm ERROR: file open error: '//trim(rpointer_fn))

        read(fid,*,iostat=ierr)tmp
        call if_true_then_msg_stop(ierr.ne.0,'norcpm ERROR: file read error: '//trim(rpointer_fn))

        close(fid)
        fn = trim(tmp)
        return
    end function current_restart_fn

    subroutine if_true_then_msg_stop(flag,msg)
        logical,intent(in) :: flag
        character(len=*), intent(in) :: msg
        if(.not.flag)return
        write(iulog,'(a)')msg
        call shr_sys_abort(msg)
    end subroutine if_true_then_msg_stop

end module norcpm_enkf
