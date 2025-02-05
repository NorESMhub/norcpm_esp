!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Tools for NorCPM On-The-Fly data assimilation
!!  Specific to gather blom data
!!  Updates:
!!      2024-08-26, Ping-Gin. Created
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module norcpm_otf
    !! blom modules
    use mod_state
    use mod_types, only: r8
    !! use captial to avoid '-Dmod_xc=mod_xc_norcpm' in buildlib
    !! so this will be the mod_xc in BLOM instead of in micom_init
    use MOD_XC, only: xcaget,xcaput, nbdy  
                 
    use dimensions, only: itdm,jtdm, idm, jdm, kdm
    use seq_comm_mct, only: seq_comm_setptrs, seq_comm_iamin  &
                        & , seq_comm_inst, seq_comm_iam &
                        & , ALLESPID, OCNID, seq_comm_mpicom
    use mpi, only: MPI_ALLGATHER, MPI_SEND, MPI_IRECV, MPI_INTEGER  &
                & ,MPI_ANY_SOURCE, MPI_REAL8, MPI_STATUS_SIZE   &
                & ,MPI_BCAST,MPI_COMM_SIZE,MPI_BARRIER, MPI_REAL4
    use netcdf, only: NF90_OPEN,NF90_INQ_VARID,NF90_GET_VAR,NF90_CLOSE &
                    ,NF90_STRERROR,NF90_NOWRITE &
                    ,NF90_INQ_DIMID,NF90_INQUIRE_DIMENSION
    !! type for p_prep_obs() and obs_readobs()
    use mod_measurement, only: measurement, OBSTYPESTRLEN
    implicit none
    save
    private 

    !! variables should be set in esp_comp_mct.F90: esp_init_mct()
    character(len=16),public :: inst_name,inst_suffix
    integer,public :: inst_index, mype ! instance index and pe(ALLESPID) of current process
    integer,public :: ninst,NINST_ESP  ! number of instances
    integer,public :: iulog            ! log file unit
    integer,public :: espcomm          ! communicator across ESPs(ALLESPID)

    !! init in init_norcpm_otf()
    integer :: nowmyiter(10) = 0
    logical,public:: amiinocn  !! init_norcpm_otf()
    integer,public :: myOCNID  !! init_norcpm_otf(), OCNID(my inst of OCN)
    integer,public :: totalpe  !! all pes of ESPs(ALLESPID)

    !! for transfering tmpX5.uf from calc_X5() to update_fields() in m_local_analysis.F90
    real(4), allocatable, public :: X5_otf(:,:,:,:) !! (nrens,nrens,ni,nj); 10x10x360x385 x4bytes ~ 53MB
    integer,public :: ufcomm,ufrank !! for update_fields(), set in EnKF.F90
    integer :: x5_m1m2(2)
    integer,allocatable :: x5_m1m2_all(:,:)

    !! for p_prep_obs() and obs_readobs(), avoid transfer data via obsvations.uf
    type(measurement), allocatable :: obs_otf(:)
    integer :: nobs_otf
    logical :: first_sync_obs_records = .true.   !! setup measurement MPI type
    integer :: mpitype_measurements      !! init when first trans
    integer :: mpitype_measurement
    logical :: debug = .false.


    !! 1st,2nd layer of dp. Use only in m_read_EN4_profile.F90
    !! only allocate it at master of ALLESP; init at norcpm_enkf.F90
    real(8), public, allocatable :: mxd1(:,:),mxd2(:,:)

    !! for perturb temp on startup run
    !! set in esp_init_mct()
    logical, public :: PERTURB_TEMP

    !! public subroutine/function
    public :: init_norcpm_otf
    public :: gather_blom_1i1l
    public :: gather_blom_1i1l_allo
    public :: put_blom_1i1l
    public :: who_need_data_raise_hand
    public :: get_pe_start_end_iter
    public :: reset_get_pe_start_end_iter
    public :: init_micom_init
    !! for calc_X5() and update_fileds()
    public :: set_uf_comm, free_uf_comm
    public :: get_x5_otf_j,x5_otf_get_m1m2,x5_otf_save_j,gather_x5_otf_to_pe0
    public :: bcast_x5_otf_from_pe0

    public :: read_nc_file_2d_and_bcast_r8,read_nc_file_3d_and_bcast_r8
    !!!! for obs_otf operation
    public :: get_obs_records_and_allocate,sync_obs_records,reset_obs_record
    public :: append_obs_record
    !!!! for do_perturb_temp in norcpm_enkf.F90
    public :: perturb_temp_with_random_seed
    !!!! to output fld data after update_fields() in EnKF.F90
    public :: save_nc_fld_update_fields,save_2d_nc_fn_vn
    !!!! test xcaput() and xcaget() in MOD_XC
    public :: test_xca_put_get
    !!!! for replace of get_climato_dim() in m_get_micom_dim.F90
    public :: read_dim3_and_bcast,read_1d_and_bcast

contains

subroutine init_norcpm_otf()
    !! init some variables
    integer :: i,ierr
    do i = 1, ninst
        amiinocn = seq_comm_iamin(OCNID(i))
        if (amiinocn)exit 
    end do ! do i = 1, ninst
    if (amiinocn) myOCNID = OCNID(i) 
    call MPI_COMM_SIZE(espcomm,totalpe,ierr) ! all esp
    if(ierr.ne.0)then
        print*,'init_norcpm_otf(): MPI_COMM_SIZE() error: ',ierr
        stop
    end if
end subroutine init_norcpm_otf


subroutine gather_blom_1i1l_allo(srcinst,varname,lev,dstrank,outdata)
    !! allocate outdata before gather_blom_1i1l()
    implicit none
    integer, intent(in)         :: srcinst  !! member to gather
    character(len=*), intent(in):: varname  !! variable gather
    integer, intent(in)         :: lev      !! level of variable
    integer, intent(in)         :: dstrank  !! return data to rank of ALLESPID
                                            !! if <0, broadcast to all
    real(r8), allocatable,  intent(out) :: outdata(:,:) !! 1 instance 1 layer 
    if (.not.allocated(outdata))allocate(outdata(itdm,jtdm))
    call gather_blom_1i1l(srcinst,varname,lev,dstrank,outdata)
end subroutine gather_blom_1i1l_allo
subroutine gather_blom_1i1l(srcinst,varname,lev,dstrank,outdata)
    !! gather BLOM data from all tasks
    !!  need be called by all cpu
    !! procedure(slow but easier for now):
    !!  1. use BLOM:xcaget() gather 1 layer data for 1 instance (myOCNID)
    !!  2. send all instance data to dstrank (ALLESPID)
    !!     or broadcast to all if dstrank .lt.0
    implicit none
    integer, intent(in)         :: srcinst  !! member to gather
    character(len=*), intent(in):: varname  !! variable gather
    integer, intent(in)         :: lev      !! level of variable
    integer, intent(in)         :: dstrank  !! return data to rank of ALLESPID
                                            !! if <0, broadcast to all
    real(r8), intent(out) :: outdata(:,:) !! 1 instance 1 layer 

    !! local
    integer  :: ocn_rank, allesp_rank
    integer  :: ocn_rootpe
    real(r8) :: tiled_var(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !! 1 tile 1 layer
    integer,allocatable  :: ocn_master(:)
    integer :: tag !! for mpi tag
    integer :: ierr, mpistat(MPI_STATUS_SIZE)
    integer ::  i
    logical :: is_srcOCN

    tag = inst_lev_var_2tag(srcinst,varname,lev)
    !allesp_rank = seq_comm_iam(ALLESPID)
    allesp_rank = mype
    ocn_rank = -1
    is_srcOCN = seq_comm_iamin(OCNID(srcinst))
    if (is_srcOCN) ocn_rank = seq_comm_iam(myOCNID)
    !! gathering 1 layer for each instance, to rank=1 of myOCNID
    !! ocn_rank .lt.0 means no ocn data in this rank
    if (is_srcOCN)then
        call varname_to_var_r8_1tile(varname,lev,tiled_var)
        if(debug) &
            print'(a,3(i0,x),a,3(x,i0))','gather_blom_1i1l(): iam, ocn_rank,inst,vn,lev, count(tiled.eq.0),size(tiled): '&
                                ,allesp_rank,ocn_rank,srcinst,varname,lev,count(tiled_var.eq.0),size(tiled_var)
        !! OCN_comm is in mod_xc, 0 is bcast, 1+ is proc (node)
        call xcaget(outdata,tiled_var,1) 

        if(debug) &
            print'(a,3(i0,x),a,3(x,i0))','gather_blom_1i1l(): iam, ocn_rank,inst,vn,lev count(outdata),size(outdata): '&
                                ,allesp_rank,ocn_rank,srcinst,varname,lev,count(outdata.eq.0),size(outdata)
    end if
    !call MPI_BARRIER(espcomm,ierr)
    !! send to dstrank of ALLESPID
    if (dstrank .lt. 0 )then !! broadcast to all pes
        !! find srcrank and bcast
        if (ocn_rank.eq.0.and.is_srcOCN)then !! only one pe
            ocn_rootpe = allesp_rank
            do i = 0,totalpe-1
                if(allesp_rank.ne.i) then
                    call MPI_SEND(ocn_rootpe,1,MPI_INTEGER,i,tag+5,espcomm,ierr)
                    if(ierr.ne.0) then
                        print*,'ERROR,gather_blom_data():115 error: ',ierr
                        stop
                    end if
                end if
            end do !! i = 0,totalpe-1
        else !! all others
            call MPI_RECV(ocn_rootpe,1,MPI_INTEGER,MPI_ANY_SOURCE,tag+5,espcomm,mpistat,ierr)
            if(ierr.ne.0) then
                print*,'ERROR,gather_blom_data():123 error: ',ierr
                stop
            end if
        end if
        if(debug) &
            print'(a,i0,x,i0)','gather_blom_1i1l(): iam, bcast from: ',allesp_rank,ocn_rootpe
        call MPI_BCAST(outdata,size(outdata),MPI_REAL8,ocn_rootpe,espcomm,ierr)
        if(debug) &
            print'(a,i0,x,i0,x,i0)','gather_blom_1i1l(): iam, num outdata.eq.0: ',allesp_rank,count(outdata.eq.0),size(outdata)
        if(ierr.ne.0) then
            print*,'ERROR,gather_blom_data(): bcast error'
            stop
        end if
        return
    end if !! broadcast to all pes then return
    !!! 1 exception, if data already at dst rank
    if(debug)print'(a,2(i0,x),a,2(i0,x),L2)','gather_blom_1i1l():184; iam, ocn_rank,vn,lev,dstrank,is_srcOCN: ' &
            ,allesp_rank,ocn_rank,varname,lev,dstrank,is_srcOCN
    if(dstrank.eq.allesp_rank .and. ocn_rank.eq.0 .and. is_srcOCN) then
        return
    end if

    !!! send outdata from ocn_rank=0,inst_index=srcinst to allesp_rank=dstrank
    ierr = 0
    if (allesp_rank .eq. dstrank)then  !! only one rank will be matched
        if(debug)print'(a,2(i0,x),a,2(i0,x),i0)','gather_blom_1i1l():193; iam, ocn_rank,vn,lev,srcinst,dstrank: ' &
            ,allesp_rank,ocn_rank,varname,lev,srcinst,dstrank
        call MPI_RECV(outdata,size(outdata),MPI_REAL8,MPI_ANY_SOURCE,tag,espcomm,mpistat,ierr)
        if(ierr.ne.0)print*,'ERROR,gather_blom_data(): recv error'
    end if
    if (ocn_rank.eq.0.and.is_srcOCN)then !! only one rank will be matched
        if(debug)print'(a,2(i0,x),a,2(i0,x),i0)','gather_blom_1i1l():198; iam, ocn_rank,vn,lev,srcinst,dstrank: ' &
            ,allesp_rank,ocn_rank,varname,lev,srcinst,dstrank
        call MPI_SEND(outdata,size(outdata),MPI_REAL8,dstrank,tag,espcomm,ierr)
        if(ierr.ne.0)print*,'ERROR,gather_blom_data(): send error'
    end if
    if(ierr.ne.0) stop
    return
end subroutine gather_blom_1i1l
subroutine put_blom_1i1l(dstinst,varname,lev,srcrank,putdata)
    !! put BLOM data from 1 to all tasks
    !!  need be called by all cpu
    !! procedure(slow but easier for now):
    !!  1. transfer 1 layer from srcrank (ALLESPID) to dstinst master
    !!  2. use BLOM:xcaput() scatter to dstinst (myOCNID)
    implicit none
    integer, intent(in)         :: dstinst  !! member to gather
    character(len=*), intent(in):: varname  !! variable gather
    integer, intent(in)         :: lev      !! level of variable
    integer, intent(in)         :: srcrank  !! return data to rank of ALLESPID
    real(r8), intent(out)       :: putdata(itdm,jtdm) !! 1 instance 1 layer 

    !! local
    integer  :: ocn_rank, allesp_rank
    integer  :: ocn_rootpe, dstocn_comm
    real(r8) :: tiled_var(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !! 1 tile 1 layer
    integer,allocatable  :: ocn_master(:)
    integer :: tag !! for mpi tag
    integer :: ierr, mpistat(MPI_STATUS_SIZE)
    logical :: is_dstOCN

    tiled_var  = -1.  !! init

    tag = inst_lev_var_2tag(dstinst,varname,lev)
    allesp_rank = mype !seq_comm_iam(ALLESPID)
    ocn_rank = -1
    is_dstOCN = seq_comm_iamin(OCNID(dstinst))
    if (is_dstOCN) ocn_rank = seq_comm_iam(OCNID(dstinst))

    !! send data from srcrank to ocn rootpe of dstinst
    ierr = 0
    if (allesp_rank .eq. srcrank)then  !! only one rank will be matched
        !! get dstinst ocn rootpe first
        !print*,'put_blom_1i1l():261, mype, send to... who? :',mype
        if (ocn_rank.eq.0.and.is_dstOCN)then !! if same pe 
            !print*,'put_blom_1i1l():263, mype, I am the one!! :',mype
        else
            call MPI_RECV(ocn_rootpe,1,MPI_INTEGER,MPI_ANY_SOURCE,tag+5,espcomm,mpistat,ierr)
            !print*,'put_blom_1i1l():263, mype, send to rank:',mype,ocn_rootpe
            call MPI_SEND(putdata,size(putdata),MPI_REAL8,ocn_rootpe,tag,espcomm,ierr)
            if(ierr.ne.0)print*,'ERROR,put_blom_1i1l(): recv error'
        end if
    else if (ocn_rank.eq.0.and.is_dstOCN)then !! only one rank will be matched
        !print*,'put_blom_1i1l():267, m pe,rneed data from:',mype,srcrank
        call MPI_SEND(allesp_rank,1,MPI_INTEGER,srcrank,tag+5,espcomm,ierr)
        !print*,'put_blom_1i1l():269, mype, recv from rank:',mype,srcrank
        call MPI_RECV(putdata,size(putdata),MPI_REAL8,MPI_ANY_SOURCE,tag,espcomm,mpistat,ierr)
        if(ierr.ne.0)print*,'ERROR,put_blom_1i1l(): send error'
    else
        putdata = -999.  !! test
    end if
    if(ierr.ne.0) stop 1

    if (is_dstOCN) then
        !! scatter to ocn tasks of dstinst
        !! mod_xc mnproc start from 1
        !print*,'put_blom_1i1l():278, before xcaput(), mype:',mype
        call xcaput(putdata,tiled_var,1)  !! it ignore some land tile
        dstocn_comm = seq_comm_mpicom(OCNID(dstinst))
        !! assign to mod_state variable
        call assign_fld_by_varname_r8_1tile(varname,lev,tiled_var)
    end if
    !print*,'put_blom_1i1l():283, mype, done:',mype
    return
end subroutine put_blom_1i1l

subroutine get_pe_start_end_iter(dstrank,myi0,myi1,maxiter,id,i0,i1,isme)
    !! bcast mi0,mi1 to all pes
    integer, intent(in)  :: dstrank    !! iter for this rank
    integer, intent(in)  :: myi0,myi1  !! iter of called rank
    integer, intent(in)  :: maxiter    !! max iteration once
    integer, intent(in)  :: id !! which id, for continue with different calls
    integer, intent(out) :: i0,i1
    logical, intent(out) :: isme
    integer :: i0i1(2),ierr
    isme = .false.
    if(dstrank.eq.mype) then !! capable continue run
        isme = .true.
        if (nowmyiter(id).eq.myi1)then !! already last iteration
            i0i1 = 0
        else
            i0i1(1) = max(nowmyiter(id)+1,myi0)
            i0i1(2) = min(i0i1(1)+maxiter-1,myi1)
            nowmyiter(id) = i0i1(2)
        end if
    end if
    call MPI_BCAST(i0i1,2,MPI_INTEGER,dstrank,espcomm,ierr)
    if(ierr.ne.0)then
        print'(a)','get_pe_start_end_iter(): MPI_BCAST() ERROR.'
        stop
    end if
    i0 = i0i1(1)
    i1 = i0i1(2)
end subroutine get_pe_start_end_iter
subroutine reset_get_pe_start_end_iter()
    nowmyiter = 0
end subroutine reset_get_pe_start_end_iter

subroutine test_xca_put_get(dstinst)
    integer, intent(in) :: dstinst
    real(r8) :: var(itdm,jtdm), var_b(itdm,jtdm)
    real(r8) :: tiled_var(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
    integer :: i,j, ocn_rank, count0, countd
    logical :: is_dstOCN

    if(.not.amiinocn) return

    ocn_rank = -1
    is_dstOCN = seq_comm_iamin(OCNID(dstinst))
    if (is_dstOCN) ocn_rank = seq_comm_iam(OCNID(dstinst))

    !! init test data, no 0 in array
    if(ocn_rank.eq.0)then
        do j = 1,jtdm
        do i = 1,itdm
            var(i,j) = (j-1)*itdm + i
        end do !i
        end do !j
        var_b = var
    end if
    !! test xcaput()
    if (is_dstOCN) then
        call xcaput(var,tiled_var,1)  !! 
        count0 = count(tiled_var.eq.0)
        if(count0.gt.0)then
            print'(a,i4,i4,i2,i7)','test_xca_put_get(): xcaput() got 0; mype,ocn_rank,dstinst,count0= ', mype,ocn_rank,dstinst,count0
        end if
    end if
    !! test xcaget()
    if (is_dstOCN) then
        call xcaget(var,tiled_var,1)  !! 
    end if
    if(ocn_rank.eq.0)then
        countd = count(var.ne.var_b)
        if(countd.gt.0)then
            print'(a,i4,i4,i2,i7)','test_xca_put_get(): xcaget() got dif; mype,ocn_rank,dstinst,countd= ', mype,ocn_rank,dstinst,countd
            call save_2d_nc_fn_vn('test_xca_put_get-0','testvar',var_b,.true.)
            call save_2d_nc_fn_vn('test_xca_put_get-1','testvar',var  ,.true.)
            stop 1
        end if
    end if

end subroutine test_xca_put_get


function who_need_data_raise_hand(raise,tag) result(dstrank)
    !! called by allesp comm, and return the rank which raise is .true.
    logical, intent(in) :: raise
    integer, intent(in) :: tag !! for mpi tag
    integer             :: dstrank
    !! local
    integer :: i,ierr,nrank, request,mpistatus
    call MPI_COMM_SIZE(espcomm, nrank, ierr)
    if(ierr.ne.0)print'(a)','ERROR,who_need_data_raise_hand(): MPI_COMM_SIZE() &
                            &failed'
    if(raise)then
        do i = 0, nrank-1
            if(i.ne.mype) then
                call MPI_ISEND(mype,1,MPI_INTEGER,i,tag,espcomm,request,ierr)
                if(ierr.ne.0)  &
                    print'(a)','ERROR,who_need_data_raise_hand(): MPI_ISEND() &
                            &failed'
            end if
        end do
        dstrank = mype
    else
        call MPI_RECV(dstrank,1,MPI_INTEGER,MPI_ANY_SOURCE,espcomm,mpistatus,ierr)
        if(ierr.ne.0) print'(a)','ERROR,who_need_data_raise_hand(): MPI_RECV() &
                    &failed'
    end if
    call MPI_BARRIER(espcomm,ierr)
    return
end function who_need_data_raise_hand

function inst_lev_var_2tag(inst,varname,lev) result(tag)
    !! combine these arguments to a integer tag for MPI
    integer,intent(in) :: inst,lev
    character(len=*),intent(in) :: varname
    integer :: tag
    !! local
    integer :: i
    logical :: uniq
    uniq = .true.
    if (uniq)then !! for a unique tag for every layer/instance/variable
        !! hope it will not overflow
        !! assume varname is only a-z
        tag = 0
        do i = 1, len(trim(varname))
            tag = tag + (ICHAR(varname(i:i))-96)*(25*(i-1))
        end do
        tag = tag*100 + inst
        tag = tag*100 + lev
    else !! just a number of sum of arguments
        tag = 0
        do i = 1, len(trim(varname))
            tag = tag + (ICHAR(varname(i:i))-96)
        end do
        tag = tag+ inst
        tag = tag+ lev
    end if
end function inst_lev_var_2tag

subroutine varname_to_var_r8_1tile(varname,lev,var)
    !! convert variable name to variable data
    !! not a good idea, but i have no better solution
    !! see NorESM/components/blom/phy/mod_state.F90
    !! AND, kfpla is integer, need be process somewhere else
    character(len=*),intent(in) :: varname
    integer, intent(in)         :: lev
    real(r8),intent(out)        :: var(:,:)
    SELECT CASE(trim(varname))
        CASE('u') !! (,,kdm*2)
            var = u(:,:,lev)
        CASE('v')
            var = v(:,:,lev)
        CASE('dp')
            var = dp(:,:,lev)
        CASE('dpu')
            var = dpu(:,:,lev)
        CASE('dpv')
            var = dpv(:,:,lev)
        CASE('temp')
            var = temp(:,:,lev)
        CASE('saln')
            var = saln(:,:,lev)
        CASE('sigma')
            var = sigma(:,:,lev)
        CASE('uflx')
            var = uflx(:,:,lev)
        CASE('vflx')
            var = vflx(:,:,lev)
        CASE('utflx')
            var = utflx(:,:,lev)
        CASE('vtflx')
            var = vtflx(:,:,lev)
        CASE('usflx')
            var = usflx(:,:,lev)
        CASE('vsflx')
            var = vsflx(:,:,lev)
        CASE('p') !! (,,kdm+1)
            var = p(:,:,lev)
        CASE('pu')
            var = pu(:,:,lev)
        CASE('pv')
            var = pv(:,:,lev)
        CASE('phi')
            var = phi(:,:,lev)
        CASE('ubflxs') !!(,,3)
            var = ubflxs(:,:,lev)
        CASE('vbflxs')
            var = vbflxs(:,:,lev)
        CASE('ub') !! (,,2)
            var = ub(:,:,lev)
        CASE('vb')
            var = vb(:,:,lev)
        CASE('pb')
            var = pb(:,:,lev)
        CASE('pbu')
            var = pbu(:,:,lev)
        CASE('pbv')
            var = pbv(:,:,lev)
        CASE('ubflxs_p')
            var = ubflxs_p(:,:,lev)
        CASE('vbflxs_p')
            var = vbflxs_p(:,:,lev)
        CASE('pb_p') !! (,,)
            var = pb_p(:,:)
        CASE('pbu_p')
            var = pbu_p(:,:)
        CASE('pbv_p')
            var = pbv_p(:,:)
        CASE('ubcors_p')
            var = ubcors_p(:,:)
        CASE('vbcors_p')
            var = vbcors_p(:,:)
        CASE('sealv')
            var = sealv(:,:)
        CASE DEFAULT
            print*,'ERROR,varname_to_var_r8_1tile(): '//trim(varname)//' not available'
            stop
    END SELECT
end subroutine varname_to_var_r8_1tile

subroutine assign_fld_by_varname_r8_1tile(varname,lev,var)
    !! assign var(:,:) to BLOM field in mod_state
    character(len=*),intent(in) :: varname
    integer, intent(in)         :: lev
    real(r8),intent(in)        :: var(:,:)

    ! local
    real(r8),allocatable :: diffld(:,:) !! difference for diagnostic

    !! allocate
    diffld = var

    SELECT CASE(trim(varname))
        CASE('u') !! belows are (,,kdm*2)
            call assign_fld_with_check(u(:,:,lev),var,diffld,lev,varname)
        CASE('v')
            call assign_fld_with_check(v(:,:,lev),var,diffld,lev,varname)
        CASE('dp')
            call assign_fld_with_check(dp(:,:,lev),var,diffld,lev,varname)
        CASE('dpu')
            call assign_fld_with_check(dpu(:,:,lev),var,diffld,lev,varname)
        CASE('dpv')
            call assign_fld_with_check(dpv(:,:,lev),var,diffld,lev,varname)
        CASE('temp')
            call assign_fld_with_check(temp(:,:,lev),var,diffld,lev,varname)
        CASE('saln')
            call assign_fld_with_check(saln(:,:,lev),var,diffld,lev,varname)
        CASE('sigma')
            call assign_fld_with_check(sigma(:,:,lev),var,diffld,lev,varname)
        CASE('uflx')
            call assign_fld_with_check(uflx(:,:,lev),var,diffld,lev,varname)
        CASE('vflx')
            call assign_fld_with_check(vflx(:,:,lev),var,diffld,lev,varname)
        CASE('utflx')
            call assign_fld_with_check(utflx(:,:,lev),var,diffld,lev,varname)
        CASE('vtflx')
            call assign_fld_with_check(vtflx(:,:,lev),var,diffld,lev,varname)
        CASE('usflx')
            call assign_fld_with_check(usflx(:,:,lev),var,diffld,lev,varname)
        CASE('vsflx')
            call assign_fld_with_check(vsflx(:,:,lev),var,diffld,lev,varname)
        CASE('p') !! belows are (,,kdm+1)
            call assign_fld_with_check(p(:,:,lev),var,diffld,lev,varname)
        CASE('pu')
            call assign_fld_with_check(pu(:,:,lev),var,diffld,lev,varname)
        CASE('pv')
            call assign_fld_with_check(pv(:,:,lev),var,diffld,lev,varname)
        CASE('phi')
            call assign_fld_with_check(phi(:,:,lev),var,diffld,lev,varname)
        CASE('ubflxs') !! belows are (,,3)
            call assign_fld_with_check(ubflxs(:,:,lev),var,diffld,lev,varname)
        CASE('vbflxs')
            call assign_fld_with_check(vbflxs(:,:,lev),var,diffld,lev,varname)
        CASE('ub') !! belows are (,,2)
            call assign_fld_with_check(ub(:,:,lev),var,diffld,lev,varname)
        CASE('vb')
            call assign_fld_with_check(vb(:,:,lev),var,diffld,lev,varname)
        CASE('pb')
            call assign_fld_with_check(pb(:,:,lev),var,diffld,lev,varname)
        CASE('pbu')
            call assign_fld_with_check(pbu(:,:,lev),var,diffld,lev,varname)
        CASE('pbv')
            call assign_fld_with_check(pbv(:,:,lev),var,diffld,lev,varname)
        CASE('ubflxs_p')
            call assign_fld_with_check(ubflxs_p(:,:,lev),var,diffld,lev,varname)
        CASE('vbflxs_p')
            call assign_fld_with_check(vbflxs_p(:,:,lev),var,diffld,lev,varname)
        CASE('pb_p') !! belows are (,,)
            call assign_fld_with_check(pb_p(:,:),var,diffld,lev,varname)
        CASE('pbu_p')
            call assign_fld_with_check(pbu_p(:,:),var,diffld,lev,varname)
        CASE('pbv_p')
            call assign_fld_with_check(pbv_p(:,:),var,diffld,lev,varname)
        CASE('ubcors_p')
            call assign_fld_with_check(ubcors_p(:,:),var,diffld,lev,varname)
        CASE('vbcors_p')
            call assign_fld_with_check(vbcors_p(:,:),var,diffld,lev,varname)
        CASE('sealv')
            call assign_fld_with_check(sealv(:,:),var,diffld,lev,varname)
        CASE DEFAULT
            print*,'ERROR,assign_fld_by_var,diffld,lev,varname)name_r8_1tile(): '//trim(varname)//' not available'
            stop
    END SELECT
    !! output diag?

    deallocate(diffld)
end subroutine assign_fld_by_varname_r8_1tile

subroutine assign_fld_with_check(dstfld, newfld, diffld, lev, vn)
    !! check and assign newfld to dstfld
    real(r8),intent(inout)  :: dstfld(:,:)
    real(r8),intent(in)     :: newfld(:,:)
    real(r8),intent(out)    :: diffld(:,:)
    integer, intent(in)     :: lev !! for dp
    character(len=*),intent(in),optional :: vn

    !! local
    logical, allocatable :: valmask(:,:)
    real(r8) :: checkv
    integer :: i,j
    logical :: check_val
    check_val = .false.
    !! test for skip dp
    !!if(trim(vn).eq.'dp')return
    !! prepare mask
    !!!! dp
    if(trim(vn).eq.'dp')then
        valmask = newfld .lt.1.e30
    else
        valmask = dp(:,:,lev) .gt.0
        valmask = valmask .and. dp(:,:,lev).lt.1.e30
    end if

    !!!! exclude values are way too large (maybe land?)
    valmask = valmask .and. dstfld.lt.1.e30
    valmask = valmask .and. newfld.lt.1.e30

    !!!! exclude values are nan
    valmask = valmask .and. .not.isnan(dstfld)
    valmask = valmask .and. .not.isnan(newfld)

    !!!! value < 0 with temp and saln
    if(trim(vn).eq.'temp' .or.trim(vn).eq.'saln')then
        valmask = valmask .and. newfld.ge.0.
    end if

    !! get difference
    diffld = 0.
    where(valmask) diffld=newfld-dstfld
    if(check_val)then
        checkv = maxval(abs(diffld))
        if(checkv.ge.10)then
            print'(a,x,i0,x,a,x,i0,x,g0)','assign_fld_with_check(): diffld is too large, mype,vn,lev,maxval = ',mype,trim(vn),lev,checkv
            do j = 1,jdm
            do i = 1,idm
                if(.not.valmask(i,j))cycle
                if(abs(diffld(i,j)).ge.10)then
                    print'(a,x,i4,x,a,x,i3,2i3,x,g0,x,g0)','    assign_fld_with_check(): mype,vn,lev,i,j,dst,new = ',mype,trim(vn),lev,i,j,dstfld(i,j),newfld(i,j)
                end if
            end do !i
            end do !j
        end if !! checkv.ge.10
    end if !! check_val

    !! assign new to dst
    where(valmask)dstfld = newfld
end subroutine assign_fld_with_check

subroutine perturb_temp_with_random_seed(seed,maxamp)
    integer, intent(in)         :: seed
    real(kind=8), intent(in)    :: maxamp
    !local
    integer, allocatable :: seed_put(:)
    integer              :: seed_size
    integer  :: i,j,k ! counters
    real(r8) :: ran(idm,jdm)
    real(r8) :: tiled_var(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

    if (.not.amiinocn) return

    !write(iulog,*)'perturb_temp_with_random_seed(), mype,seed = ',mype,seed
    !temp = temp + float(mype)*0.001
    !return

    call random_seed(size=seed_size)
    allocate(seed_put(seed_size))

    do k = 1,kdm
        seed_put = seed*100+k
        call random_seed(put=seed_put)
        call random_number(ran(1,1)) !! init
        call varname_to_var_r8_1tile('temp',k,tiled_var)
        call random_number(ran) !! max is 1.0
        !! for positive and negative
        !write(*,'(a,i0,x,5f7.4,4x,5f7.4)')'perturb_temp_with_random_seed(): seed, ran(1:5,1) = ',seed,ran(1:5,1), ((ran(1:5,1)-0.5)*maxamp*2)
        tiled_var(1:idm,1:jdm) = tiled_var(1:idm,1:jdm) + ((ran-0.5)*maxamp*2)
        call assign_fld_by_varname_r8_1tile('temp',k,tiled_var)
    end do ! k

end subroutine perturb_temp_with_random_seed

subroutine init_micom_init(imem,imemmax,enssize)
    integer, intent(out) :: imem,imemmax,enssize
    imem = seq_comm_inst(myOCNID)
    enssize = ninst
end subroutine init_micom_init

subroutine set_uf_comm(color)
    !! set a temp comm for update_field()
    integer, intent(in) :: color
    integer :: ierr
    call MPI_COMM_SPLIT(espcomm,color,mype,ufcomm,ierr)
    if(ierr.ne.0)print*,'set_uf_comm(): MPI_COMM_SPLIT(): ',ierr
    call MPI_COMM_RANK(ufcomm,ufrank,ierr)
    if(ierr.ne.0)print*,'set_uf_comm(): MPI_COMM_RANK(): ',ierr
end subroutine set_uf_comm
subroutine free_uf_comm()
    integer :: ierr
    call MPI_COMM_FREE(ufcomm,ierr)
    ufrank = -1
end subroutine free_uf_comm
subroutine x5_otf_get_m1m2(m1,m2,ni,nj)
    !! save full X5_otf for now
    integer, intent(in) :: m1,m2,ni,nj
    integer :: ierr
    x5_m1m2(1) = m1
    x5_m1m2(2) = m2
    if(allocated(X5_otf))deallocate(X5_otf)
    allocate(X5_otf(ninst,ninst,ni,nj))
    if(allocated(x5_m1m2_all))deallocate(x5_m1m2_all)
    allocate(x5_m1m2_all(2,totalpe))
    !! gather m1m2
    !print*,'x5_otf_get_m1m2():710; mype = ',mype
    call MPI_ALLGATHER(x5_m1m2,2,MPI_INTEGER, &
        x5_m1m2_all,totalpe*2,MPI_INTEGER,espcomm,ierr)
end subroutine x5_otf_get_m1m2
subroutine x5_otf_save_j(x5,j)
    real(4),intent(in) :: x5(:,:,:)
    integer,intent(in) :: j
    if(mype.eq.0)then
        X5_otf(:,:,:,j) = x5
    else
        X5_otf(:,:,:,j-x5_m1m2(1)+1) = x5
    end if
end subroutine x5_otf_save_j
subroutine gather_x5_otf_to_pe0()
    !! gather X5_otf and bcast to all
    integer :: i,j,tag,n
    integer :: ierr, mpistat(MPI_STATUS_SIZE)
    n = size(X5_otf(:,:,:,1))
    !! gathers, maybe change to allgatherv?
    if(mype.eq.0)then
        do i = 2, totalpe
            if(x5_m1m2_all(1,i).gt.x5_m1m2_all(2,i))cycle
            do j = x5_m1m2_all(1,i),x5_m1m2_all(2,i)
                tag = i-1+j
                !print*,'gather_x5_otf_to_pe0():732; j,n = ',j,n
                call MPI_RECV(X5_otf(:,:,:,j),n,MPI_REAL4,i-1,tag,espcomm,mpistat,ierr)
            end do ! j
        end do !i
    else
        if(x5_m1m2(1).le.x5_m1m2(2))then
            do j = x5_m1m2(1),x5_m1m2(2)
                tag = mype+j
                call MPI_SEND(X5_otf(:,:,:,j-x5_m1m2(1)+1),n,MPI_REAL4,0,tag,espcomm,ierr)
            end do !j
        end if
    end if
end subroutine gather_x5_otf_to_pe0
subroutine bcast_x5_otf_from_pe0()
    integer :: ierr
    call MPI_BCAST(X5_otf,size(X5_otf),MPI_REAL4,0,espcomm,ierr)  
    !print*,'bcast_x5_otf_from_pe0(): mype,X5_otf.eq.0 count = ',mype,count(X5_otf.eq.0)
end subroutine bcast_x5_otf_from_pe0
subroutine get_x5_otf_j(j,x5)
    integer, intent(in) :: j
    real(4),intent(out) :: x5(:,:,:)
    x5 = X5_otf(:,:,:,j)
end subroutine get_x5_otf_j


subroutine read_nc_file_2d_and_bcast_r8(fn,vn,var)
    !! read real8 data variable in filename
    !! then use mpi_bcast() to allesp comm.
    !! should be called by all proc in espcomm
    !! not done yet
    character(len=*),intent(in)  :: fn,vn
    real(kind=8),intent(out)     :: var(:,:)
    ! local
    integer :: varid, ncid, ierr
    !! read from netcdf file (assume netcdf is not compiled with parallel)
    !! also assume the dimension is matched
    if(mype.eq.0)then
        call donc(NF90_OPEN(trim(fn),NF90_NOWRITE,ncid))
        call donc(NF90_INQ_VARID(ncid,trim(vn),varid))
        call donc(NF90_GET_VAR(ncid,varid,var))
        call donc(NF90_CLOSE(ncid))
    end if
    !! bcast to members in allesp comm
    call MPI_BCAST(var,size(var),MPI_REAL8,0,espcomm,ierr)
    if(ierr.ne.0)then
        write(iulog,'(a,i0)')'read_nc_file_2d_and_bcast_r8() error: ',ierr
        stop ierr
    end if
end subroutine read_nc_file_2d_and_bcast_r8
subroutine read_nc_file_3d_and_bcast_r8(fn,vn,var)
    !! read real8 data variable in filename
    !! then use mpi_bcast() to allesp comm.
    !! should be called by all proc in espcomm
    !! not done yet
    character(len=*),intent(in)  :: fn,vn
    real(kind=8),intent(out)     :: var(:,:,:)
    ! local
    integer :: varid, ncid, ierr
    !! read from netcdf file (assume netcdf is not compiled with parallel)
    !! also assume the dimension is matched
    if(mype.eq.0)then
        call donc(NF90_OPEN(trim(fn),NF90_NOWRITE,ncid))
        call donc(NF90_INQ_VARID(ncid,trim(vn),varid))
        call donc(NF90_GET_VAR(ncid,varid,var))
        call donc(NF90_CLOSE(ncid))
    end if
    !! bcast to members in allesp comm
    call MPI_BCAST(var,size(var),MPI_REAL8,0,espcomm,ierr)
    if(ierr.ne.0)then
        write(iulog,'(a,i0)')'read_nc_file_3d_and_bcast_r8() error: ',ierr
        stop ierr
    end if
end subroutine read_nc_file_3d_and_bcast_r8
subroutine donc(ierr)
    integer,intent(in) :: ierr
    !! handle error from netcdf
    if(ierr.ne.0) then
        write(iulog,'(2a)')'netcdf error: ',NF90_STRERROR(ierr)
        stop 1
    end if
end subroutine donc


subroutine reset_obs_record()
    !! clean up the obs data
    if(allocated(obs_otf))deallocate(obs_otf)
end subroutine reset_obs_record

subroutine append_obs_record(newobs)
    !! replace the write_wet_file() in p_prep_obs.F90
    type(measurement),intent(in) :: newobs(:)
    !! local 
    type(measurement),allocatable :: tmp(:)
    integer :: n_otf, n_new
    if(allocated(obs_otf))then
        n_otf = size(obs_otf)
        n_new = size(newobs)
        tmp = obs_otf
        deallocate(obs_otf)
        allocate(obs_otf(n_otf+n_new))
        obs_otf(:n_otf) = tmp
        obs_otf(n_otf+1:) = newobs
        deallocate(tmp)
    else
        obs_otf = newobs
    end if
end subroutine append_obs_record

subroutine get_obs_records_and_allocate(obs,nobs)
    !! replace the reading observation.uf in EnKF.F90
    type(measurement),optional,allocatable,intent(out) :: obs(:)
    integer, optional, intent(out) :: nobs
    if(present(obs))then
        if(allocated(obs))deallocate(obs)
        obs = obs_otf
    end if
    if(present(nobs))nobs = size(obs_otf)
end subroutine get_obs_records_and_allocate

subroutine sync_obs_records()
    !! gather and sync all read obs
    integer :: provide_nobs(totalpe), nobs, i, ier, n1
    type(measurement),allocatable :: newobs(:)

    if (first_sync_obs_records)then !! setup MPI type
        call set_mpi_measurement()  !! output is mpitype_measurements
        first_sync_obs_records = .false.
    end if
    !! gather nobs from all proc
    nobs = 0
    if(allocated(obs_otf))nobs = size(obs_otf)
    call MPI_ALLGATHER(nobs,1,MPI_INTEGER, &
        provide_nobs,totalpe,MPI_INTEGER,espcomm,ier)
    !! bcast new obs
    allocate(newobs(sum(provide_nobs)))
    if(nobs.gt.0)then
        n1 = 1
        if(mype.gt.0) n1 = sum(provide_nobs(:mype))+1
        newobs(n1:n1+nobs-1) = obs_otf
    end if
    !!!! it should be compliable to parallel reading, 
    !!!! but norcpm_enkf.F90 is sequency now.
    n1 = 1
    if(.true.)then !! for parallel reading, not ready
        do i = 1, totalpe
            if(provide_nobs(i).eq.0)cycle
            call MPI_BCAST(newobs(n1:n1+provide_nobs(i)-1),provide_nobs(i) &
                            ,mpitype_measurements,i-1,espcomm,ier)
            n1 = n1+provide_nobs(i)
        end do ! i = 1, totalpe
    end if
    !! assign newobs to nobs_otf
    call reset_obs_record()
    obs_otf = newobs
    deallocate(newobs)

end subroutine sync_obs_records

subroutine set_mpi_measurement()
    !! setup mpi type for measurement
    !! output is mpitype_measurements
    use mpi, only: MPI_GET_ADDRESS,MPI_ADDRESS_KIND, &
        MPI_TYPE_CREATE_STRUCT,MPI_TYPE_CREATE_RESIZED,MPI_TYPE_COMMIT, &
        MPI_CHARACTER,MPI_LOGICAL
    type(measurement) :: dum(2)
    integer, parameter :: natt = 19
    integer(kind=MPI_ADDRESS_KIND) :: offsets(natt),extent
    integer :: atypes(natt), alens(natt)
    integer :: ier, i
    character,allocatable :: msg
    !write(*,'(a,i0,a)',advance='no')'set_mpi_measurement(), pe=',mype,'... '
    !! get attrs address related to var
    msg = 'MPI_GET_ADDRESS():674 error'
        call MPI_GET_ADDRESS(dum(1)%d          ,offsets( 1), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%var        ,offsets( 2), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%id         ,offsets( 3), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%lon        ,offsets( 4), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%lat        ,offsets( 5), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%depth      ,offsets( 6), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%ipiv       ,offsets( 7), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%jpiv       ,offsets( 8), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%ns         ,offsets( 9), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%a1         ,offsets(10), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%a2         ,offsets(11), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%a3         ,offsets(12), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%a4         ,offsets(13), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%status     ,offsets(14), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%i_orig_grid,offsets(15), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%j_orig_grid,offsets(16), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%h          ,offsets(17), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%date       ,offsets(18), ier);call err_mpi(ier,msg)
        call MPI_GET_ADDRESS(dum(1)%orig_id    ,offsets(19), ier);call err_mpi(ier,msg)
    offsets(2:) = offsets(2:) - offsets(1)
    offsets(1) = 0
    !! set type and length
        i= 1; atypes(i) = MPI_REAL8;    alens(i) = 1  ! d
        i= 2; atypes(i) = MPI_REAL8;    alens(i) = 1  ! var
        i= 3; atypes(i) = MPI_CHARACTER;alens(i) = OBSTYPESTRLEN  ! id
        i= 4; atypes(i) = MPI_REAL8;    alens(i) = 1  ! lon
        i= 5; atypes(i) = MPI_REAL8;    alens(i) = 1  ! lat
        i= 6; atypes(i) = MPI_REAL8;    alens(i) = 1  ! depth
        i= 7; atypes(i) = MPI_INTEGER;  alens(i) = 1  ! ipiv
        i= 8; atypes(i) = MPI_INTEGER;  alens(i) = 1  ! jpiv
        i= 9; atypes(i) = MPI_INTEGER;  alens(i) = 1  ! ns
        i=10; atypes(i) = MPI_REAL8;    alens(i) = 1  ! a1
        i=11; atypes(i) = MPI_REAL8;    alens(i) = 1  ! a2
        i=12; atypes(i) = MPI_REAL8;    alens(i) = 1  ! a3
        i=13; atypes(i) = MPI_REAL8;    alens(i) = 1  ! a4
        i=14; atypes(i) = MPI_LOGICAL;  alens(i) = 1  ! status
        i=15; atypes(i) = MPI_INTEGER;  alens(i) = 1  ! i_orig_grid
        i=16; atypes(i) = MPI_INTEGER;  alens(i) = 1  ! j_orig_grid
        i=17; atypes(i) = MPI_REAL8;    alens(i) = 1  ! h
        i=18; atypes(i) = MPI_INTEGER;  alens(i) = 1  ! date
        i=19; atypes(i) = MPI_INTEGER;  alens(i) = 1  ! orig_id
    !! create MPI type for single measurement
    call MPI_TYPE_CREATE_STRUCT(natt,alens,offsets,atypes,mpitype_measurement,ier)
    call err_mpi(ier,'MPI_TYPE_CREATE_STRUCT():717')

    !! get len between elements in arrary. offsets is reused
    call MPI_GET_ADDRESS(dum(1)%d ,offsets(1), ier)
    call err_mpi(ier,'MPI_GET_ADDRESS():721')
    call MPI_GET_ADDRESS(dum(2)%d ,offsets(2), ier)
    call err_mpi(ier,'MPI_GET_ADDRESS():723')
    extent = offsets(2) - offsets(1)
    !! create resized, not sure what it's means
    call MPI_TYPE_CREATE_RESIZED(mpitype_measurement &
        ,0_MPI_ADDRESS_KIND,extent,mpitype_measurements,ier)
    call err_mpi(ier,'MPI_TYPE_CREATE_RESIZED():728')
    call MPI_TYPE_COMMIT(mpitype_measurements,ier)
    call err_mpi(ier,'MPI_TYPE_COMMIT():730')
    !write(*,*)'done.',mpitype_measurements
end subroutine set_mpi_measurement

subroutine err_mpi(ierr,msg)
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: msg
    character(len=255) :: errstr
    integer :: n,i
    if (ierr.ne.0) then
        call MPI_ERROR_STRING(ierr,errstr,n,i)
        print*,msg,' ',errstr(1:n)
        stop
    end if
end subroutine err_mpi

subroutine save_2d_nc_fn_vn(fn,vn,var,addpe)
    use netcdf, only: nf90_noclobber,nf90_create,nf90_def_dim,nf90_real8, &
                      nf90_def_var,nf90_enddef,nf90_put_var,nf90_eexist
    real(8),intent(in) :: var(:,:)
    character(len=*),intent(in) :: fn,vn
    logical, intent(in) :: addpe

    character(len=255) :: ncfn
    integer :: ierr,ncid,dims(2),xdimid,ydimid,varid

    dims = shape(var)
    if(addpe)then
        write(ncfn,'(2a,i4.4,a)')fn,'-pe',mype,'.nc'
    else
        ncfn = fn
    end if

    ierr = nf90_create(trim(ncfn),nf90_noclobber,ncid)
    if(ierr.eq. nf90_eexist)then
        print*,'save_2d_nc_fn_vn(): file exists, skip: ',trim(ncfn)
        return
    end if
    ierr = nf90_def_dim(ncid,'x',dims(1),xdimid)
    ierr = nf90_def_dim(ncid,'y',dims(2),ydimid)
    ierr = nf90_def_var(ncid,vn,nf90_real8,[xdimid,ydimid],varid)
    ierr = nf90_enddef(ncid)
    ierr = nf90_put_var(ncid,varid,var)
    ierr = nf90_close(ncid)
    write(*,'(2a)')'save_2d_nc_fn_vn(): ',trim(ncfn)
end subroutine save_2d_nc_fn_vn

subroutine save_nc_fld_update_fields(fld,vn,lev,tag)
    real(8),intent(in) :: fld(:,:)
    character(len=*),intent(in) :: vn,tag
    integer, intent(in) :: lev
    !! local
    character(len=100) :: ncfn

    write(ncfn,'(2a,i4.4,3a,i2.2,a)')trim(tag),'-',mype,'-',trim(vn),'-',lev,'.nc'
    call save_2d_nc_fn_vn(ncfn,vn,fld,.false.)
end subroutine save_nc_fld_update_fields

subroutine read_dim3_and_bcast(fn,nx,ny,nz)
    !! replacment of get_climato_dim() in m_get_micom_dim.F90
    character(len=*),intent(in) :: fn
    integer,intent(out) :: nx,ny,nz
    integer :: mpibuffer(3),ierr
    integer :: ncid,dimid
    if(mype.eq.0)then
        call donc(NF90_OPEN(trim(fn),NF90_NOWRITE,ncid))
        call donc(NF90_INQ_DIMID(ncid,'x',dimid))
        call donc(NF90_INQUIRE_DIMENSION(ncid,dimid,len=nx))
        call donc(NF90_INQ_DIMID(ncid,'y',dimid))
        call donc(NF90_INQUIRE_DIMENSION(ncid,dimid,len=ny))
        call donc(NF90_INQ_DIMID(ncid,'depth',dimid))
        call donc(NF90_INQUIRE_DIMENSION(ncid,dimid,len=nz))
        call donc(NF90_CLOSE(ncid))
    end if
    mpibuffer = [nx,ny,nz]
    call MPI_BCAST(mpibuffer,3,MPI_INTEGER,0,espcomm,ierr)
    nx = mpibuffer(1)
    ny = mpibuffer(2)
    nz = mpibuffer(3)
end subroutine read_dim3_and_bcast
subroutine read_1d_and_bcast(fn,vn,n,var)
    !! replacment of get_climato_dim() in m_get_micom_dim.F90
    character(len=*),intent(in) :: fn,vn
    integer,intent(in) :: n
    real(kind=8),intent(out) :: var(n)
    integer :: ierr
    integer :: ncid,varid
    if(mype.eq.0)then
        call donc(NF90_OPEN(trim(fn),NF90_NOWRITE,ncid))
        call donc(NF90_INQ_VARID(ncid,vn,varid))
        call donc(NF90_GET_VAR(ncid,varid,var))
        call donc(NF90_CLOSE(ncid))
    end if
    call MPI_BCAST(var,n,MPI_REAL8,0,espcomm,ierr)
end subroutine read_1d_and_bcast


end module norcpm_otf
