#!/bin/bash -e

case=noresm2_betzy_25_5inst  ## 4 inst, use modified compset
    NINST=5
    NINST_ESP=5
    compset=NHISTNCPM
    NTASKS=$((128*$NINST))   ## 8 nodes for 4 inst with MULTI_DRIVER=FALSE
    NTASKS_OCN=$((123*$NINST)) ## 123 * 4
        #### blom_dimensions: Available processor counts: 32 42 63 77 91 123 156 186 256 354
    NTASKS_ESP=$((128*$NINST_ESP))
    MULTI_DRIVER=FALSE
        ## TRUE:  total tasks = NTASKS*NINST,
        ## FALSE: total tasks = NTASKS, share to NINST

pecount=S
compset=${compset:-NHIST} ## NorESM
res=${res:-f19_tn14}  ## NorESM
noresmdir="/cluster/projects/nn9039k/people/pgchiu/NorESM_git/NorESM"
caseroot=.

## cleaning old build
echo "To clean old build:"
echo "    rm -rf ~/work/{archive,noresm,noresm2_cases}/${case}"
if [ "$1" == "delete" ] ;then
    echo 'cleaning old bld'
    rm -rf ~/work/{archive,noresm,noresm2_cases}/${case}
    rm -rf ${caseroot}/${case}/
    if [ "$2" == "only" ] ;then
        exit
    fi
fi

${noresmdir}/cime/scripts/./create_newcase --case ${case} --mach betzy --res ${res} --compset ${compset} --pecount $pecount --project nn9039k


cd ${caseroot}/${case}

if [ ! -z "$NTASKS" ] ; then
    ./xmlchange NTASKS_ATM=$NTASKS
    ./xmlchange ROOTPE_ATM=0
    
    ./xmlchange NTASKS_LND=$NTASKS
    ./xmlchange ROOTPE_LND=0
    
    ./xmlchange NTASKS_ICE=$NTASKS
    ./xmlchange ROOTPE_ICE=0
    
    ./xmlchange NTASKS_OCN=$NTASKS_OCN
    #### blom_dimensions: Available processor counts: 32 42 63 77 91 123 156 186 256 354
    ./xmlchange ROOTPE_OCN=0
    
    ./xmlchange NTASKS_CPL=$NTASKS
    ./xmlchange ROOTPE_CPL=0
    
    ./xmlchange NTASKS_GLC=$NTASKS
    ./xmlchange ROOTPE_GLC=0
    
    ./xmlchange NTASKS_ROF=$NTASKS
    ./xmlchange ROOTPE_ROF=0
    
    ./xmlchange NTASKS_WAV=$NTASKS
    ./xmlchange ROOTPE_WAV=0

    ./xmlchange NTASKS_ESP=$NTASKS_ESP
    ./xmlchange ROOTPE_ESP=0
fi

## pause
./xmlchange NINST_ATM=$NINST
./xmlchange NINST_LND=$NINST
./xmlchange NINST_ICE=$NINST
./xmlchange NINST_ROF=$NINST
./xmlchange NINST_GLC=$NINST
./xmlchange NINST_OCN=$NINST
./xmlchange NINST_WAV=$NINST
./xmlchange NINST_ESP=$NINST_ESP
test "$NINST" -gt 1 && test ! -z "$MULTI_DRIVER" && ./xmlchange MULTI_DRIVER=$MULTI_DRIVER ## or out of memory
## pause ocn only, per day
test "$NINST" -gt 1 && ./xmlchange PAUSE_ACTIVE_OCN=TRUE,PAUSE_ACTIVE_ICE=TRUE,PAUSE_N=1,PAUSE_OPTION=nday,ESP_RUN_ON_PAUSE=TRUE

./xmlchange STOP_N=5
./xmlchange STOP_OPTION=nday
./xmlchange JOB_WALLCLOCK_TIME=01:00:00
./xmlchange RUN_STARTDATE=1970-01-01

## soft link makes life easier
ln -s $(./xmlquery RUNDIR --value)      ./run
ln -s $(./xmlquery EXEROOT --value)     ./bld
ln -s $(./xmlquery DOUT_S_ROOT --value) ./archive

./case.setup
#./preview_namelists
./case.build
./case.submit
