!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module
!
! Structure:
!
! Revisions:
!   2023-11 Ping-Gin Chiu: File created
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module norcpm_enkf

    use m_EnKF, only: EnKF
    use m_prep_obs, only: p_prep_obs
    use m_micom_ensemble_init, only: micom_ensemble_init
    use m_ensave, only: p_ensave
    use mpi
    use netcdf
    use shr_sys_mod, only : shr_sys_abort
    use ESMF, only: ESMF_Clock
    use seq_timemgr_mod, only: seq_timemgr_EClockGetData
    use perf_mod, only: t_startf,t_stopf !!,t_initf,t_finalizef

    implicit none
    save
    private
    logical :: first = .true.
    integer :: iulog
    integer :: mpicomm, mype, ninst
    integer :: mpisize
    integer :: rootpe = 0
    logical :: master = .false.

    !! namelist
    character(len=255)             :: INPUTDATA,OCNGRIDFILE,MEAN_MOD_DIR,RES
    character(len=255),dimension(3) :: PRODUCERLIST, OBSLIST, REF_PERIODLIST, MONTHLY_ANOM, &
                                       COMBINE_ASSIM
    integer :: ANOMALYASSIM, ANOM_CPL, OSAS, RFACTOR, ENSAVE, ENSSIZE
    character(len=4)   :: fforano
    character(len=255) :: ANALYSIS_DIRNAME,RESULT_DIRNAME
    !! namelist modelio
    character(len=255) :: diri,diro,logfile

    public :: norcpm_assim_step

    !! the order of variable in namelist is meanful, need be same as esp_in
    namelist /norcpm/INPUTDATA,OCNGRIDFILE &
                    ,PRODUCERLIST, OBSLIST, REF_PERIODLIST, MONTHLY_ANOM, COMBINE_ASSIM &
                    ,ANOMALYASSIM, ANOM_CPL, OSAS, RFACTOR,fforano, ANALYSIS_DIRNAME,RESULT_DIRNAME &
                    ,RES,MEAN_MOD_DIR,ENSAVE,ENSSIZE
    namelist /modelio/diri,diro,logfile

contains

    subroutine norcpm_assim_step(mpicomm,EClock)
        integer, intent(in) :: mpicomm !! ESP MPI comm
        type(ESMF_Clock),intent(in) :: EClock
        logical :: should_return = .false.
        integer :: ierr, EnKF_CNT = 0,fid
        character(len=2) :: ec_str !! str for EnKF_CNT
        integer :: year,month,day,sec
        integer :: i
        integer :: commsize

        if (first) call init_save_variables()
        call MPI_BARRIER(mpicomm,ierr)

        call seq_timemgr_EClockGetData(EClock, year,month,day,sec)

        should_return = .false.
        if (master)then
            if(is_restarts_EnKF_already(current_restart_fn('ocn',1))) then
                write(iulog,*)'norcpm: Restart files are already EnKF, skip.'
                should_return = .true.
            end if
        end if
        call MPI_BCAST(should_return,1,MPI_LOGICAL,rootpe,mpicomm,ierr)
        if(should_return)return

        call link_forecasts()
        call MPI_BARRIER(mpicomm,ierr)

        !! prep_obs*3, distribute tasks, will it work?
        !! need more examine assim_step.sh and EnKF.F90
        !! for now it's plain imply
        !! need review
        EnKF_CNT = 0
        do i = 1, size(OBSLIST) 
            call prep_prep_obs(trim(OBSLIST(i)),trim(PRODUCERLIST(i)) &
                ,trim(REF_PERIODLIST(i)) ,trim(COMBINE_ASSIM(i)),year,month)
            call MPI_BARRIER(mpicomm,ierr)
            if(master)then
                call exec_in_analysis_dir('prep_obs')
                call exec_in_analysis_dir( 'mv observations.uf observations.uf_' & 
                            //trim(OBSLIST(i))//'.'//trim(PRODUCERLIST(i)))
            end if

            call MPI_BARRIER(mpicomm,ierr)
            if (trim(COMBINE_ASSIM(i)).ne.'1')cycle

            EnKF_CNT = EnKF_CNT+1
            if(master)then
                write(iulog,'(a,i2)')'norcpm: EnKF_CNT= ',EnKF_CNT
                flush(iulog)
                write(ec_str,'(i2)')EnKF_CNT
                ec_str = adjustl(ec_str)
                call write_analysisfields_file(trim(ANALYSIS_DIRNAME)//'/analysisfields.in')
                call write_enkf_prm_file(ENSSIZE,RFACTOR,trim(ANALYSIS_DIRNAME)//'/enkf.prm')
                call exec_in_analysis_dir( &
                    & 'cat observations.uf_* > observations.uf ; &
                    &  rm -f observations.uf_*')
            end if
            if(ENSAVE.eq.1)then
                call exec_in_analysis_dir('ensave','forecast')
                if(master .and. ENSAVE.eq.1)then
                    call exec_in_analysis_dir('mv forecast_avg.nc forecast_avg_'//trim(ec_str)//'.nc')
                    ! no aicen in BLOM, skip ensave_ice()
                end if
            end if
            call MPI_BARRIER(mpicomm,ierr)
            call exec_in_analysis_dir('EnKF')
            call MPI_BARRIER(mpicomm,ierr)
            if(master)then
                call exec_in_analysis_dir( '&
                    & mv enkf_diag.nc enkf_diag_'//trim(ec_str)//'.nc; &
                    & mv tmpX5.uf tmpX5_'//trim(ec_str)//'.uf ')
            end if
            if(ENSAVE.eq.1)then
                call exec_in_analysis_dir('ensave','forecast')
                if(master) &
                    call exec_in_analysis_dir('mv forecast_avg.nc analysis_avg_'//trim(ec_str)//'.nc')
                ! no aicen in BLOM, skip ensave_ice()
            end if
            if(master)write(iulog,'(a)')'norcpm: Finished with EnKF; call number :'//trim(ec_str)
        end do

        ! archive assimilation files, need revise
        if(master)then
            !!call exec_in_analysis_dir('mkdir -p '//RESULT_DIRNAME//'/')
            call archive_assimilation(year,month)
        end if

        call MPI_BARRIER(mpicomm,ierr)
        call exec_in_analysis_dir('micom_ensemble_init')
        call MPI_BARRIER(mpicomm,ierr)

        !call EXECUTE_COMMAND_LINE
            !('cd '//ANALYSIS_DIRNAME//'; &
            !& rm -f forecast???.nc forecast_ice???.nc aiceold???.nc viceold???.nc &
            !&       observations.uf enkf.prm* infile.data*')
    end subroutine norcpm_assim_step

    subroutine archive_assimilation(year,month)
        integer, intent(in) :: year,month
        character(len=4) :: yyyy
        character(len=2) :: mm
        write(yyyy,'(i4.4)')year
        write(mm,  '(i2.2)')month
        call exec_in_analysis_dir('mkdir -p '//trim(RESULT_DIRNAME)//'/'//yyyy//'_'//mm)
        call exec_in_analysis_dir('mv enkf_diag_*.nc observations-*.nc tmpX5_*.uf '//trim(RESULT_DIRNAME)//'/'//yyyy//'_'//mm)
        if(ENSAVE.eq.1)then
            call exec_in_analysis_dir('mv analysis_*avg_*.nc forecast_*avg_*.nc '//trim(RESULT_DIRNAME)//'/'//yyyy//'_'//mm)
            call exec_in_analysis_dir('mv SAL.nc TEM.nc '//trim(RESULT_DIRNAME)//'/'//yyyy//'_'//mm)
        end if
        !! no aicen and vicen in BLOM restart, need revise
    end subroutine archive_assimilation

    subroutine exec_in_analysis_dir(cmd,arg1,arg2)
        character(len=*),intent(in) :: cmd
        character(len=*),intent(in),optional :: arg1,arg2
        character(:), allocatable :: precmd
        character(len=255) :: rundir
        integer :: ierr

        if (master)write(iulog,'(2a)')'norcpm exec_in_analysis_dir(): ',cmd
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

        if(trim(cmd).eq.'micom_ensemble_init')then
            call t_startf('micom_ensemble_init')
            call chdir(trim(ANALYSIS_DIRNAME))
            call micom_ensemble_init(ENSSIZE)
            call chdir(rundir)
            call t_stopf('micom_ensemble_init')
            return
        end if

        precmd = 'cd '//trim(ANALYSIS_DIRNAME)//'/ ;'
        call EXECUTE_COMMAND_LINE(precmd//' '//cmd, exitstat=ierr)
        if(ierr.ne.0) call shr_sys_abort('norcpm ERROR: cmd: '//precmd//' '//cmd)
    end subroutine exec_in_analysis_dir

    subroutine prep_prep_obs(OBSTYPE,PRODUCER,REF_PERIOD,COMB_ASSIM,year,month,day,sec)
        character(len=*),intent(in) :: OBSTYPE,PRODUCER,REF_PERIOD,COMB_ASSIM
        integer, intent(in)         :: year,month
        integer, intent(in),optional:: day,sec
        character(len=255) :: fn,fn2,dst,msg,cmd
        character(len=2) :: mm
        character(len=4) :: yyyy
        write(mm,'(i2.2)')month
        write(yyyy,'(i4.4)')year

        if(mype.eq.0) then ! obs data
            fn = trim(INPUTDATA)//"/Obs/"//OBSTYPE//'/'//PRODUCER//'/'//yyyy//'_'//mm//'.nc'
            fn2= trim(INPUTDATA)//"/Obs/"//OBSTYPE//'/'//PRODUCER//'/'//yyyy//'_'//mm//'_pre.nc'
            msg= 'norcpm ERROR: Obs data missing: '//trim(fn)
            dst= trim(ANALYSIS_DIRNAME)//'/'//yyyy//'_'//mm//'.nc'
            call link_if_exist_or_stop(fn=trim(fn),fn2=trim(fn2),dst=trim(dst),stopmsg=trim(msg))
        end if

        if(mype.eq.1) then ! model mean
            fn = trim(MEAN_MOD_DIR)//'/Free-average'//mm//'-'//REF_PERIOD//'.nc'
            msg= 'norcpm ERROR: model mean data missing: '//trim(fn)
            dst= trim(ANALYSIS_DIRNAME)//'/mean_mod.nc'
            call link_if_exist_or_stop(fn=trim(fn),dst=trim(dst),stopmsg=trim(msg))
        end if

        if(mype.eq.2) then ! obs unc
            fn = trim(INPUTDATA)//'/Input/NorESM/'//trim(RES)//'/'//PRODUCER//'/'//&
                &trim(RES)//'_'//OBSTYPE//'_obs_unc_'//fforano//'.nc'
            dst= trim(ANALYSIS_DIRNAME)//'/obs_unc_'//OBSTYPE//'.nc'
            msg= 'norcpm: obs unc missing, but maybe fine: '//trim(fn)
            call link_if_exist_or_stop(fn=trim(fn),dst=trim(dst),stopmsg=trim(msg),stoprun=.false.)
        end if

        if(mype.eq.3) then ! obs mean
            fn = trim(INPUTDATA)//'/Obs/'//OBSTYPE//'/'//PRODUCER//'/'//OBSTYPE//'_avg_'//mm//'-'//REF_PERIOD//'.nc'
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
            fn = trim(INPUTDATA)//'/Input/EnKF/infile.data.'//trim(OBSTYPE)//'.'//trim(PRODUCER)
            dst = trim(ANALYSIS_DIRNAME)//'/infile.data'
            cmd = 'sed -e "s/yyyy/'//yyyy//'/" -e "s/mm/'//mm//'/" "'//trim(fn)//'" > "'//trim(dst)//'"'
            write(iulog,'(a)')'norcpm: run cmd='//trim(cmd)
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
        write(iulog,'(a)')'norcpm: soft link: '//fn//' -> '//dst
        if(is_file_exist(fn))then
            write(iulog,'(a)')"ln -sfr "//fn//' '//dst
            call EXECUTE_COMMAND_LINE("ln -sfr "//fn//' '//dst,exitstat=ierr)
            if (ierr.eq.0)return
        end if
        if(present(fn2))then
            if(is_file_exist(fn2))then
                write(iulog,'(a)')"ln -sfr "//fn2//' '//dst
                call EXECUTE_COMMAND_LINE("ln -sfr "//fn2//' '//dst,exitstat=ierr)
                if (ierr.eq.0)return
            end if
        end if
        if (ierr.ne.0)then
            call get_environment_variable(name='PWD',value=pwd,status=ierr)
            write(iulog,'(2a)')'norcpm ERROR: link process error,pwd=',trim(pwd)
            write(iulog,'(2a)')'norcpm ERROR: link process error,src=',trim(fn)
            write(iulog,'(2a)')'norcpm ERROR: link process error,dst=',trim(dst)
        end if
        if(stoprun_.and.present(stopmsg)) call shr_sys_abort(stopmsg)
        if(stoprun_) call  shr_sys_abort('norcpm ERROR: link_if_exist_or_stop()')
    end subroutine link_if_exist_or_stop

    logical function is_file_exist(fn)
        character(len=*),intent(in) :: fn
        inquire(file=fn,exist=is_file_exist)
    end function is_file_exist

    subroutine link_forecasts()
        integer :: i
        character(len=3) :: i3
        if(mype .le. (ninst-1))then
            write(i3,'(i3.3)')mype+1
            call delink_and_link(current_restart_fn('ocn',mype+1) &
                ,trim(ANALYSIS_DIRNAME)//'/forecast'//i3//'.nc')
            return
        end if
        if(.false.)then !! skip ice for now
        if(mype .le. (ninst*2-1))then
            write(i3,'(i3.3)')mype-ninst+1
            call delink_and_link(current_restart_fn('ice',mype-ninst+1) &
                ,trim(ANALYSIS_DIRNAME)//'/forecast_ice'//i3//'.nc')
            if(.false.)then !! no aicen and vicen in BLOM restart, need revise
                call extract_ncfile(trim(ANALYSIS_DIRNAME)//'/forecast_ice'//i3//'.nc' &
                    ,trim(ANALYSIS_DIRNAME)//'/aiceold'//i3//'.nc','aicen')
                call extract_ncfile(trim(ANALYSIS_DIRNAME)//'/forecast_ice'//i3//'.nc' &
                    ,trim(ANALYSIS_DIRNAME)//'/viceold'//i3//'.nc','vicen')
            end if
        end if
        end if
    end subroutine link_forecasts

    subroutine delink_and_link(src,dst)
        character(len=*), intent(in) :: src,dst
        character(len=255) :: cmd
        integer :: exitstat
        cmd = "test $(stat -c '%h' '"//src//"') -eq 1 && test ! -h '"//src//"'"
        write(iulog,'(a)')trim(cmd)
        call EXECUTE_COMMAND_LINE(cmd,exitstat=exitstat)
        if(exitstat.ne.0)then !! if src reference to other file
            call EXECUTE_COMMAND_LINE('cp "'//src//'" "'//src//'_dst"')
            call EXECUTE_COMMAND_LINE('mv --backup=t "'//src//'" "'//src//'_orig"')
            call EXECUTE_COMMAND_LINE('mv "'//src//'_dst" "'//src//'"')
        end if

        call EXECUTE_COMMAND_LINE("ln -sfr "//src//' '//dst)
    end subroutine delink_and_link

    subroutine init_save_variables()
        integer :: ierr

        call MPI_COMM_RANK(mpicomm,mype,ierr)
        call MPI_COMM_SIZE(mpicomm,mpisize,ierr)
        call if_true_then_msg_stop(ierr.ne.MPI_SUCCESS, &
                    'norcpm ERROR: MPI_COMM_RANK() error')
        if(mype .eq. rootpe)master = .true.

        call read_namelist_norcpm('esp_in')

        if(master) call EXECUTE_COMMAND_LINE('rm -rf   '//trim(ANALYSIS_DIRNAME))
        if(master) call EXECUTE_COMMAND_LINE('mkdir -p '//trim(ANALYSIS_DIRNAME))

        write(iulog,'(a,i4,i4)')'norcpm: init, mype,mpisize: ',mype,mpisize
        first = .false.
    end subroutine init_save_variables

    subroutine read_namelist_norcpm(fn)
        character(len=*) :: fn
        integer :: fid, ios
        if (.true.)then !! the namelist generate is not ready yet, check buildnml
            open(file=fn,action='read',newunit=fid,iostat=ios)
            if(ios.ne.0) write(iulog,*)"norcpm, Error: "//trim(fn)//' namelist open error'
            read(nml=norcpm,iostat=ios,unit=fid)
            if(ios.ne.0) write(iulog,*)"norcpm, Error: "//trim(fn)//' namelist read error'
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
        write(iulog,*)'norcpm: write_enkf_prm_file(): ENSSIZE = ',ENSSIZE
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
        write(fid,'(a)')'u         1 53'
        write(fid,'(a)')'v         1 53'
        write(fid,'(a)')'dp        1 53'
        write(fid,'(a)')'temp      1 53'
        write(fid,'(a)')'saln      1 53'
        !!write(fid,'(a)')'uflx      1 53'
        !!write(fid,'(a)')'vflx      1 53'
        !!write(fid,'(a)')'utflx     1 53'
        !!write(fid,'(a)')'vtflx     1 53'
        !!write(fid,'(a)')'usflx     1 53'
        !!write(fid,'(a)')'vsflx     1 53'
        write(fid,'(a)')'pb        1 1'
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

    logical function is_restarts_EnKF_already(ncfn)
        character(len=*), intent(in) :: ncfn
        is_restarts_EnKF_already = .false.
        write(iulog,*)'norcpm is restart file EnKFed? ',ncfn
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
            write(iulog,*)'norcpm ERROR: '//trim(NF90_STRERROR(ierr))
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
