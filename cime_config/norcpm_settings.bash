#!/bin/bash
# Experiment default setting, can be set below to override default.


## case setting
: ${RES:=CASERES}  ## resolution, f09_tn14 (MM) or f19_tn14 (LM)
: ${ENSSIZE:=NINST_ESP}
: ${MPIRUN:=mpirun}
: ${NTASKS_OCN:=NTASKS_OCN}

## DA setting
## 0: disable, 1: enable
### wait for be deleted ##: ${EnKF_Version:=1}
: ${ASSIMULATEMONTHDAY:=15}  ## day of DA each month
: ${INPUTDATA:=NORCPM_INPUTDATA}
: ${OCNGRIDFILE:=OCNGRIDFILE}
: PRODUCERLIST=${PRODUCERLIST:='HADISST2 EN421 EN421'}  ## norcpm_ana_f09_tn14 1980-2010
: OBSLIST=${OBSLIST:='SST TEM SAL'}    ## data to assimilate
: REF_PERIODLIST=${REF_PERIODLIST:='1980-2010 1980-2010 1980-2010'} ## ref. period of climatology
: MONTHLY_ANOM=${MONTHLY_ANOM:='1 1 1'}      ## Use anomaly data to 
: COMBINE_ASSIM=${COMBINE_ASSIM:='0 0 1'}      ## Run assimilation individually
: ${ANOMALYASSIM:=1}  # Anomaly assimilation, need consistant with MONTHLY_ANOM
: ${ANOM_CPL:='0'}    # anomaly coupled (Koseki et al. 2017)
: ${OSAS:='0'}        # Param estimation

## Arrays for pbs_enkf, array element number should be consistant with num of '1' in $COMBINE_ASSIM
## Variable to DA in model. Format: <varname> <1 or 0> <levels>
ANALYSIS_FIELDS='u         1 53
    v         1 53
    dp        1 53
    temp      1 53
    saln      1 53
    uflx      1 53
    vflx      1 53
    utflx     1 53
    vtflx     1 53
    usflx     1 53
    vsflx     1 53
    pb        1 1
    ub        1 1
    vb        1 1
    ubflx     1 1
    vbflx     1 1
    ubflxs    1 1
    vbflxs    1 1
    ubcors_p  0 0
    vbcors_p  0 0
    phi       0 0
    sealv     0 0
    ustar     0 0
    buoyfl    0 0 '

## Namelist for EnKF.F90
: ${RFACTOR:=1}   # Slow assimilation start
ENKF_PRM="
    &method
        methodtag = 'DEnKF'
    /
    &ensemble
        enssize = ${ENSSIZE}
    /
    &localisation
        locfuntag = 'Gaspari-Cohn'
        locrad = 1500.0
    /
    &moderation
        infl = 1.00
        rfactor1 = ${RFACTOR}
        rfactor2 = 4.0
        kfactor = 2.0
    /
    &files
    /
    &prmest /"


if  (( ${ANOMALYASSIM} ))  ; then
    fforano=anom
else
    fforano=ff
fi
#FOLLOWING are used to developing
: ${ANALYSIS_DIRNAME:=NORCPM_ANALYSIS_DIRNAME}  ## DA work dirname, located at RUNDIR/
: ${RESULT_DIRNAME:=NORCPM_RESULT_DIRNAME}      ## DA diag data dirname, located at RUNDIR/
: ${MPIRUN_ARGS:=''}
: ${NTASKS_ENKF:=128}
: ${NTASKS_MICOM_INIT:=NTASKS_MICOM_INIT}
