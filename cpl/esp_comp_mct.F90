module esp_comp_mct

    use mct_mod,          only: mct_aVect
    use esmf,             only: ESMF_Clock
    use seq_timemgr_mod,  only: seq_timemgr_EClockGetData
    use seq_cdata_mod,    only: seq_cdata,seq_cdata_setptrs
    use seq_infodata_mod, only: seq_infodata_PutData,seq_infodata_GetData, &
        seq_infodata_start_type_cont, seq_infodata_start_type_brnch, &
        seq_infodata_start_type_start
    use shr_file_mod,     only: shr_file_getLogUnit
    use shr_sys_mod,      only: shr_sys_flush, shr_sys_abort
    use norcpm_enkf,      only: norcpm_assim_step_otf, do_perturb_temp
    use norcpm_otf,       only: inst_index,mype,ninst,NINST_ESP,iulog,espcomm &
                              &,inst_name,inst_suffix, init_norcpm_otf &
                              &,PERTURB_TEMP

    implicit none
    save
    private 
    namelist /case/ NINST_ESP,PERTURB_TEMP

    !--------------------------------------------------------------------------
    ! Public interfaces
    !--------------------------------------------------------------------------

    public :: esp_init_mct
    public :: esp_run_mct
    public :: esp_final_mct
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
    subroutine esp_init_mct( EClock, cdata, x2d, d2x, NLFilename )
        use seq_comm_mct, only: seq_comm_suffix, seq_comm_inst, seq_comm_name, seq_comm_iam &
                                , num_inst_esp, ALLESPID, seq_comm_mpicom
        use seq_infodata_mod, only: seq_infodata_type
        !use spmdMod, only: masterproc, spmd_init

        ! !INPUT/OUTPUT PARAMETERS:

        type(ESMF_Clock)            , intent(inout) :: EClock
        type(seq_cdata)             , intent(inout) :: cdata
        type(mct_aVect)             , intent(inout) :: x2d, d2x
        character(len=*), optional  , intent(in)    :: NLFilename

        ! Local variables
        ! cdata
        integer :: ESPID
        integer :: mpicom
        type(seq_infodata_type), pointer :: infodata ! NorESM driver info data
        integer :: ierr
        character(len=32) :: starttype
        !EOP
        !-------------------------------------------------------------------------------


        call seq_infodata_PutData(cdata%infodata, esp_present=.true., &
             esp_prognostic=.true., esp_phase=1)
        call shr_file_getLogUnit(iulog)

        !! init instance name, index, suffix
        !! init io log unit
        !! 

        call seq_cdata_setptrs(cdata,ID=ESPID,mpicom=mpicom,infodata=infodata)

        inst_name = seq_comm_name(ESPID)
        inst_index = seq_comm_inst(ESPID)
        inst_suffix = seq_comm_suffix(ESPID)
        mype = seq_comm_iam(ALLESPID) !! only this instance

        !write(iulog,'(a,a,i4,2a,i4)') &
        !    'norcpm, ESP init, inst name,index,suffix,iam: ' &
        !    ,trim(inst_name),inst_index,' ',trim(inst_suffix),mype

        call read_namelist('esp_in')
        ninst = NINST_ESP

        espcomm = seq_comm_mpicom(ALLESPID)

        call seq_infodata_GetData(infodata, start_type=starttype)
        
        mype = seq_comm_iam(ALLESPID) !! cross instance myrank

        if(mype.eq.0)print'(2a)','esp_init_mct(), start_type = ',starttype
        if(trim(starttype).eq.trim(seq_infodata_start_type_start))then
            if(PERTURB_TEMP) do_perturb_temp = .true.
        end if
        call init_norcpm_otf() !! init something
    end subroutine esp_init_mct


    subroutine esp_run_mct( EClock, cdata, x2d, d2x)
        implicit none

        ! !INPUT/OUTPUT PARAMETERS:

        type(ESMF_Clock)            ,intent(inout) :: EClock
        type(seq_cdata)             ,intent(inout) :: cdata
        type(mct_aVect)             ,intent(inout) :: x2d
        type(mct_aVect)             ,intent(inout) :: d2x

        !EOP
        !-------------------------------------------------------------------------------

        call norcpm_assim_step_otf(espcomm,EClock)

    end subroutine esp_run_mct

    subroutine esp_final_mct( EClock, cdata, x2d, d2x)
        implicit none

        !----- arguments -----
        type(ESMF_Clock)            ,intent(inout) :: EClock
        type(seq_cdata)             ,intent(inout) :: cdata
        type(mct_aVect)             ,intent(inout) :: x2d
        type(mct_aVect)             ,intent(inout) :: d2x

        !EOP
        !-------------------------------------------------------------------------------

        print*,'norcpm, ESP final'
    end subroutine esp_final_mct

!===============================================================================

    subroutine read_namelist(fn)
        character(len=*) :: fn
        ! local
        integer :: fid, ios

        open(file=fn,action='read',newunit=fid,iostat=ios)
        read(nml=case,iostat=ios,unit=fid)
        if(ios.ne.0)then
            write(iulog,*)"norcpm, Error: "//trim(fn)//' namelist read error'
        end if
        close(fid)
    end subroutine read_namelist

end module esp_comp_mct
