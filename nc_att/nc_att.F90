!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! program: nc_att
! desc: query or set global attribute in a netcdf file
! usage:
!        nc_att <netcdf file> <attribute> [value]
!    if value is exist, set the attribute value to netcdf file
!    if value is empty, print the attribute value.
!
! reversion:
!    Aug2023, Ping-Gin: Create, to mark the restart file after EnKF
!                               and avoid NCO as dependency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program nc_att
    use netcdf
    implicit none
    character(len=300) :: ncfn, attr, val
    integer :: ncid, attnum, ios, xtype
    integer :: ival
    real    :: rval
    logical :: isexist
    logical :: debug=.true.

    call getarg(1,ncfn)
    call getarg(2,attr)
    call getarg(3,val)

    inquire(file=ncfn,exist=isexist)
    if(.not.isexist)then
        !write(*,'(2a)')'File not found, exit: ',trim(ncfn)
        stop 'nc_att ERROR: File not found, exit'
    end if

    if (trim(attr) .eq. '')then
        !write(*,'(a)')'attribute not exist, exit..'
        stop 'nc_att ERROR: attribute not exist, exit'
    end if

    if (trim(val) .eq. '')then  !! query attr in ncfn
        call check(nf90_open(ncfn,NF90_NOWRITE,ncid))
        ios  = nf90_inquire_attribute(ncid,NF90_GLOBAL,trim(attr),attnum=attnum,xtype=xtype)
        if (ios .ne. nf90_noerr) stop 
        select case (xtype)  !! should be interface, but it is a small tool...
            case (NF90_FLOAT)
                call check(nf90_get_att(ncid,NF90_GLOBAL,trim(attr),rval))
                write(*,'(f)')rval
            case (NF90_INT)
                call check(nf90_get_att(ncid,NF90_GLOBAL,trim(attr),ival))
                write(*,'(i)')ival
            case (NF90_CHAR)
                call check(nf90_get_att(ncid,NF90_GLOBAL,trim(attr),val))
                write(*,'(a)')trim(val)
            case DEFAULT
                stop 'nc_att ERROR, not suitable xtype'
        end select

    else                        !! set attr in ncfn
        call check(nf90_open(ncfn,NF90_WRITE,ncid))
        call check(nf90_redef(ncid))
        call check(nf90_put_att(ncid,NF90_GLOBAL,trim(attr),trim(val)))
    end if

    call check(nf90_close(ncid))
    
    contains
      subroutine check(status)
        integer, intent ( in) :: status
        
        if(status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop 'nc_att ERROR.'
        end if
      end subroutine check  
end program nc_att

