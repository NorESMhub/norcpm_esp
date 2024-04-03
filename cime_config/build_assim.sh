#!/bin/bash -e
    dobuild=0
    for i in EnKF ensave  ensave_ice  fixenkf_cice  micom_init  nc_att prep_obs ; do
        test -f "$ANALYSIS/$i" || dobuild=1
    done
    test "$dobuild" == 0 && exit || true

## build EnKF and others in Ingo's setup
    echo + build EnKF
    export MACH=betzy ## for default EnKF env setting
    mkdir -p $ANALYSISROOT/bld/EnKF/TMP
    cd $ANALYSISROOT/bld/EnKF
    cp -f $ASSIMROOT/shared/* . 
    cp -f $ASSIMROOT/EnKF/* . 
    make clean
    make $MAKE_J

    echo + build prep_obs
    mkdir -p $ANALYSISROOT/bld/prep_obs/TMP
    cd $ANALYSISROOT/bld/prep_obs
    cp -f $ASSIMROOT/shared/* . 
    cp -f $ASSIMROOT/prep_obs/* . 
    make clean
    make $MAKE_J

    echo + build ensave and fixenkf
    mkdir -p $ANALYSISROOT/bld/ensave_fixenkf/TMP
    cd $ANALYSISROOT/bld/ensave_fixenkf
    cp -f $ASSIMROOT/shared/* . 
    cp -f $ASSIMROOT/ensave_fixenkf/* . 
    make clean
    make $MAKE_J

    echo + build micom_init
    mkdir -p $ANALYSISROOT/bld/micom_init
    cd $ANALYSISROOT/bld/micom_init
    cp -f $ASSIMROOT/micom_init/* . 

    ## use current grid dimension 
    ##cp -f $EXEROOT/ocn/src/dimensions.F .
    ## use pre-set grid dimension (bad idea, need fix)
    cp -f $ASSIMROOT/../cime_config/dimensions.F.n32_k53_tnx1v4 ./dimensions.F

    ## it's tripolar grid
    sed -i 's/^CPPDEFS =.*/CPPDEFS = -DMPI -DARCTIC/' Makefile
    make clean
    make $MAKE_J

    echo + build tools/nc_att
    mkdir -p $ANALYSISROOT/bld/nc_att
    cd $ANALYSISROOT/bld/nc_att
    cp -f $SRCROOT/nc_att/* . 
    make clean
    make
    cp -f nc_att ../..

    cd ../..
    mkdir -p $ANALYSIS
    cp EnKF ensave  ensave_ice  fixenkf_cice  micom_init  nc_att prep_obs  $ANALYSIS/

    echo + link ocn grid.nc
    ## get grid.nc path fron blom namelist
    gridpath=$(grep 'GRFILE *=' $RUNDIR/ocn_in* |head -n1 | cut -d'"' -f2 | cut -d"'" -f2)
    ln -sf "$gridpath" $ANALYSIS/grid.nc

