#!/bin/bash -e

case=norcpm_esp_test01
    NINST=10
    compset=NHISTNCPM
    STOP_N=1
    STOP_OPTION=nmonth
    NTASKS_CPL=$((128*$NINST))   
    NTASKS_ATM=$((128*$NINST))   
    NTASKS_LND=$((128*$NINST))   
    NTASKS_ICE=$((128*$NINST))   
    NTASKS_OCN=$((123*$NINST)) #### blom_dimensions: Available processor counts: 32 42 63 77 91 123 156 186 256 354
    NTASKS_ROF=$((128*$NINST))   
    NTASKS_GLC=$((128*$NINST))   
    NTASKS_WAV=$((128*$NINST))   
    NTASKS_ESP=$((128*$NINST))
    ROOTPE_CPL=0
    ROOTPE_ATM=0
    ROOTPE_LND=0
    ROOTPE_ICE=0
    ROOTPE_OCN=0
    ROOTPE_ROF=0
    ROOTPE_GLC=0
    ROOTPE_WAV=0
    ROOTPE_ESP=0

    MULTI_DRIVER=FALSE
    ## DO NOT set anything about PAUSE

pecount=S
compset=${compset:-NHIST} ## NorESM
res=${res:-f19_tn14}  ## NorESM
noresmdir="./NorESM"
caseroot=.

## cleaning old build
echo "To clean old build:"
echo "    rm -rf ~/work/{archive,noresm,.}/${case}"
if [ "$1" == "delete" ] ;then
    echo 'cleaning old bld'
    rm -rf ~/work/{archive,noresm,.}/${case}
    rm -rf ${caseroot}/${case}/
    if [ "$2" == "only" ] ;then
        exit
    fi
fi

${noresmdir}/cime/scripts/./create_newcase --case ${case} --mach betzy --res ${res} --compset ${compset} --project nn9039k 


cd ${caseroot}/${case}

## set NTASKS and ROOTPE
    ./xmlchange NTASKS_ATM=$NTASKS_ATM,ROOTPE_ATM=$ROOTPE_ATM
    ./xmlchange NTASKS_LND=$NTASKS_LND,ROOTPE_LND=$ROOTPE_LND
    ./xmlchange NTASKS_ICE=$NTASKS_ICE,ROOTPE_ICE=$ROOTPE_ICE
    ./xmlchange NTASKS_OCN=$NTASKS_OCN,ROOTPE_OCN=$ROOTPE_OCN
    ./xmlchange NTASKS_CPL=$NTASKS_CPL,ROOTPE_CPL=$ROOTPE_CPL
    ./xmlchange NTASKS_GLC=$NTASKS_GLC,ROOTPE_GLC=$ROOTPE_GLC
    ./xmlchange NTASKS_ROF=$NTASKS_ROF,ROOTPE_ROF=$ROOTPE_ROF
    ./xmlchange NTASKS_WAV=$NTASKS_WAV,ROOTPE_WAV=$ROOTPE_WAV
    ./xmlchange NTASKS_ESP=$NTASKS_ESP,ROOTPE_ESP=$ROOTPE_ESP

## multi-instance
    ./xmlchange NINST_ATM=$NINST
    ./xmlchange NINST_LND=$NINST
    ./xmlchange NINST_ICE=$NINST
    ./xmlchange NINST_ROF=$NINST
    ./xmlchange NINST_GLC=$NINST
    ./xmlchange NINST_OCN=$NINST
    ./xmlchange NINST_WAV=$NINST
    ./xmlchange NINST_ESP=$NINST

./xmlchange STOP_N=STOP_N
./xmlchange STOP_OPTION=STOP_OPTION
./xmlchange  --subgroup case.run  JOB_WALLCLOCK_TIME=72:00:00
./xmlchange RUN_STARTDATE=1982-01-01

## soft link makes life easier
ln -s $(./xmlquery RUNDIR --value)      ./run
ln -s $(./xmlquery EXEROOT --value)     ./bld
ln -s $(./xmlquery DOUT_S_ROOT --value) ./archive

./case.setup
#./preview_namelists
./case.build
./case.submit
