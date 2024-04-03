#!/bin/sh -ev

echo + PERFORM ASSIMILATION UPDATE `date`
cd $ANALYSISROOT

#every second month we reduced rfactor by 1
nmonth=$(( (10#$yr - 10#$START_YEAR) * 12 + 10#$mm - 10#$START_MONTH )) 
RFACTOR=$(( $RFACTOR_START - $nmonth / 2 ))
if [[ $RFACTOR -lt 1 || $STOP_N_FORECAST ]]
then 
  RFACTOR=1
fi  
echo ++ SET RFACTOR TO $RFACTOR

echo ++ LINK FORECASTS
for MEMBER in `seq -w $MEMBER1 $MEMBERN`
do
  CASE=${EXPERIMENT}_$SDATE_PREFIX${SDATE}_$MEMBER_PREFIX$MEMBER
  RELPATH=$EXPERIMENT/${EXPERIMENT}_$SDATE_PREFIX$SDATE/$CASE
  EXEROOT=$WORK/noresm/$RELPATH
  if [[ $((10#$yr)) -eq $((10#$START_YEAR)) && $((10#$mm)) -eq $((10#$START_MONTH)) ]]
  then 
    if [ $REF_SUFFIX_MEMBER1 ]
    then
      REF_MEMBER1=`echo $REF_SUFFIX_MEMBER1 | tail -3c`
      if [ $RUN_TYPE == branch ]
      then
        REF_MEMBER=`printf %02d $((10#$MEMBER+10#$REF_MEMBER1-10#$MEMBER1))`
      else
        REF_MEMBER=$REF_MEMBER1
      fi
      REF_SUFFIX=`basename $REF_SUFFIX_MEMBER1 $REF_MEMBER1`$REF_MEMBER
    fi
    RESCASE=$REF_EXPERIMENT$REF_SUFFIX
    # avoid modifying restart files from reference experiment
    mv $EXEROOT/run/${RESCASE}.micom.r.$yr-$mm-15-00000.nc $EXEROOT/run/${RESCASE}.micom.r.$yr-$mm-15-00000.ncORIG
    mv $EXEROOT/run/${RESCASE}.cice.r.$yr-$mm-15-00000.nc $EXEROOT/run/${RESCASE}.cice.r.$yr-$mm-15-00000.ncORIG 
    nccopy -6 $EXEROOT/run/${RESCASE}.micom.r.$yr-$mm-15-00000.ncORIG $EXEROOT/run/${RESCASE}.micom.r.$yr-$mm-15-00000.nc
    nccopy -6 $EXEROOT/run/${RESCASE}.cice.r.$yr-$mm-15-00000.ncORIG $EXEROOT/run/${RESCASE}.cice.r.$yr-$mm-15-00000.nc
  else
    RESCASE=$CASE
    if [[ `ncdump -k $EXEROOT/run/${RESCASE}.micom.r.$yr-$mm-15-00000.nc | cut -d" " -f1` == netCDF-4 ]]
    then 
      mv $EXEROOT/run/${RESCASE}.micom.r.$yr-$mm-15-00000.nc $EXEROOT/run/${RESCASE}.micom.r.$yr-$mm-15-00000.ncORIG
      nccopy -6 $EXEROOT/run/${RESCASE}.micom.r.$yr-$mm-15-00000.ncORIG $EXEROOT/run/${RESCASE}.micom.r.$yr-$mm-15-00000.nc
    fi 
    if [[ `ncdump -k $EXEROOT/run/${RESCASE}.cice.r.$yr-$mm-15-00000.nc | cut -d" " -f1` == netCDF-4 ]]
    then 
      mv $EXEROOT/run/${RESCASE}.cice.r.$yr-$mm-15-00000.nc $EXEROOT/run/${RESCASE}.cice.r.$yr-$mm-15-00000.ncORIG 
      nccopy -6 $EXEROOT/run/${RESCASE}.cice.r.$yr-$mm-15-00000.ncORIG $EXEROOT/run/${RESCASE}.cice.r.$yr-$mm-15-00000.nc
    fi 
  fi 
  ln -sf $EXEROOT/run/${RESCASE}.micom.r.$yr-$mm-15-00000.nc forecast0${MEMBER}.nc
  ln -sf $EXEROOT/run/${RESCASE}.cice.r.$yr-$mm-15-00000.nc forecast_ice0${MEMBER}.nc
  ncks -O -v aicen forecast_ice0${MEMBER}.nc aiceold0${MEMBER}.nc
  ncks -O -v vicen forecast_ice0${MEMBER}.nc viceold0${MEMBER}.nc
done

echo ++ SKIP ASSIMILATION IF ASSIMILATION UPDATED ALREADY APPLIED TO RESTART FILES
if [[ `ncdump -h forecast0${MEMBER1}.nc | grep micom_init | wc -l` -gt 0 ]]
then
  return 0
fi 

ENKF_CNT=0 # Counter of EnKF sequential call
echo ++ PREPARE OBSERVATIONS AND DO SEQUENTIAL/CONCURRENT ASSIMILATION
OBSLIST=(${OBSLIST[*]}) # convert list to array 
PRODUCERLIST=(${PRODUCERLIST[*]})
REF_PERIODLIST=(${REF_PERIODLIST[*]})
COMBINE_ASSIM=(${COMBINE_ASSIM[*]})
for iobs in ${!OBSLIST[*]}
do
  OBSTYPE=${OBSLIST[$iobs]}
  PRODUCER=${PRODUCERLIST[$iobs]}
  REF_PERIOD=${REF_PERIODLIST[$iobs]}
  COMB_ASSIM=${COMBINE_ASSIM[$iobs]}    #sequential/joint observation assim 
  if [ -e $INPUTDATA/obs/$OBSTYPE/$PRODUCER/${yr}_${mm}.nc ]
  then  
    ln -sf $INPUTDATA/obs/$OBSTYPE/$PRODUCER/${yr}_${mm}.nc .
  elif [ -e $INPUTDATA/obs/$OBSTYPE/$PRODUCER/${yr}_${mm}_pre.nc ]
  then 
    ln -sf $INPUTDATA/obs/$OBSTYPE/$PRODUCER/${yr}_${mm}_pre.nc ${yr}_${mm}.nc
  else
    echo "$INPUTDATA/obs/$OBSTYPE/$PRODUCER/${yr}_${mm}.nc missing, we quit" ; exit 1
  fi
  if [ -e $MEAN_MOD_DIR/Free-average$mm-${REF_PERIOD}.nc ]
  then 
    ln -sf $MEAN_MOD_DIR/Free-average$mm-${REF_PERIOD}.nc mean_mod.nc 
  else
    echo "$MEAN_MOD_DIR/Free-average$mm-${REF_PERIOD}.nc missing, we quit" ; exit 1
  fi
  if [ -f $INPUTDATA/enkf/$RES/$PRODUCER/${RES}_${OBSTYPE}_obs_unc_anom.nc ]
  then
    ln -sf $INPUTDATA/enkf/$RES/$PRODUCER/${RES}_${OBSTYPE}_obs_unc_anom.nc  obs_unc_${OBSTYPE}.nc
  fi
  ln -sf $INPUTDATA/obs/$OBSTYPE/$PRODUCER/${OBSTYPE}_avg_${mm}-${REF_PERIOD}.nc mean_obs.nc || { echo "Error $INPUTDATA/obs/$OBSTYPE/$PRODUCER/${OBSTYPE}_avg_$mm-${REF_PERIOD}.nc missing, we quit" ; exit 1 ; }
  ln -sf $OCNGRIDFILE grid.nc 
  cat $INPUTDATA/enkf/infile.data.${OBSTYPE}.$PRODUCER | sed  "s/yyyy/$yr/" | sed  "s/mm/$mm/" > infile.data
  time srun -n 1 ./prep_obs
  mv observations.uf observations.uf_${OBSTYPE}.$PRODUCER
  if (( $COMB_ASSIM ))
  then
    let ENKF_CNT=ENKF_CNT+1
    cat observations.uf_* > observations.uf
    rm -f observations.uf_*
    cp -f $ASSIMROOT/analysisfields_${ENKF_CNT}.in analysisfields.in
    cat $ASSIMROOT/enkf.prm_${ENKF_CNT} | sed  "s/XXX/$RFACTOR/" > enkf.prm
    sed -i s/"enssize =".*/"enssize = "$ENSSIZE/g enkf.prm

    if (( $ENSAVE ))
    then 
      echo +++ COMPUTE PRE-ASSIMILATION ENSEMBLE MEANS FOR OCEAN AND SEA ICE 
      time ./ensave forecast $ENSSIZE 
      mv forecast_avg.nc forecast_avg_${ENKF_CNT}.nc
      if [ `grep aicen analysisfields.in | wc -l` -gt 0 ]
      then 
        time ./ensave_ice forecast_ice $ENSSIZE 
        mv forecast_ice_avg.nc forecast_ice_avg_${ENKF_CNT}.nc
      fi 
    fi

    echo +++ CALL ENKF
    time srun -n $ENKF_NTASKS ./EnKF enkf.prm
    mv enkf_diag.nc enkf_diag_${ENKF_CNT}.nc      
    mv tmpX5.uf tmpX5_${ENKF_CNT}.uf

    if (( $ENSAVE ))
    then
      echo +++ COMPUTE POST-ASSIMILATION ENSEMBLE MEANS FOR OCEAN AND SEA ICE 
      time srun -n $ENSSIZE ./ensave forecast $ENSSIZE 
      mv forecast_avg.nc analysis_avg_${ENKF_CNT}.nc
      if [ `grep aicen analysisfields.in | wc -l` -gt 0 ]
      then
        time ./ensave_ice forecast_ice $ENSSIZE 
        mv forecast_ice_avg.nc analysis_ice_avg_${ENKF_CNT}.nc
      fi
    fi

    echo 'Finished with EnKF; call number :' $ENKF_CNT
    date
  fi
done #OBS list

echo ++ ARCHIVE ASSIMILATION FILES
mkdir -p $WORK/noresm/$EXPERIMENT/${EXPERIMENT}_$SDATE_PREFIX$SDATE/RESULT/${yr}_$mm
mv enkf_diag_*.nc observations-*.nc tmpX5_*.uf $WORK/noresm/$EXPERIMENT/${EXPERIMENT}_$SDATE_PREFIX$SDATE/RESULT/${yr}_$mm
if (( $ENSAVE ))
then
  mv analysis_*avg_*.nc forecast_*avg_*.nc $WORK/noresm/$EXPERIMENT/${EXPERIMENT}_$SDATE_PREFIX$SDATE/RESULT/${yr}_$mm
fi 

if [ `grep aicen analysisfields.in | wc -l` -gt 0 ]
then
  echo ++ FIXENKF_ICE - POST-ASSIMILATION CORRECTION OF SEA ICE STATE
  for MEMBER in `seq -w $MEMBER1 $MEMBERN`
  do
    time ./fixenkf_cice $MEMBER & # process members in parallel on first node
  done
  wait # wait until all members are finished
fi 

echo ++ MICOM_INIT - POST-ASSIMILATION CORRECTION OF OCEAN STATE    
time srun -n $(( MICOM_INIT_NTASKS_PER_MEMBER * ENSSIZE )) ./micom_init $ENSSIZE

if (( $ENSAVE ))
then
  echo ++ COMPUTE FINAL ENSEMBLE MEANS FOR OCEAN AND SEA ICE 
  time srun -n $ENSSIZE ./ensave forecast $ENSSIZE
  mv forecast_avg.nc $WORK/noresm/$EXPERIMENT/${EXPERIMENT}_$SDATE_PREFIX$SDATE/RESULT/${yr}_$mm/fix_analysis_avg.nc
  if [ `grep aicen analysisfields.in | wc -l` -gt 0 ]
  then
    time ./ensave_ice forecast_ice $ENSSIZE
    mv forecast_ice_avg.nc $WORK/noresm/$EXPERIMENT/${EXPERIMENT}_$SDATE_PREFIX$SDATE/RESULT/${yr}_$mm/fix_analysis_ice_avg.nc
  fi
fi

echo ++ FINISH ASSIM POST-PROCESSING `date`
rm -f forecast???.nc forecast_ice???.nc aiceold???.nc viceold???.nc
rm -f observations.uf enkf.prm* infile.data*

echo + FINISHED ASSIMILATION UPDATE
