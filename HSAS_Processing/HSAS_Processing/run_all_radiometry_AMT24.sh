#!/bin/bash

# Hayley Evers-King 08/02/17
# Code to run all processing for AMT27 radiometry (hypersas calibrated by Tartu)

# Preprocessing required:
# Hypersas data must be calibrated outside this script if you are using a Satlantic/Satcon formatted calibration

# modified by gdal on 16/11/2018 to process amt24 data



# USER PARAMETERS
#-------------------------------------------------------------------------------

# os parameters
# date -> gdate for osx
OSX=False

# date parameters 
FIRST_DATE='2014-09-25' # this is jday=268 which is when we started running the IOP uway system
LAST_DATE='2014-11-01'
# LAST_DATE='2014-09-29'

# instrument parameters
RUN_HSAS=True

# Previously used versions to deal with different calibration versions for AMT26 - AMT 27 just using single calibration
VERSION='v1'

# L0 merging removed as only needed to AMT26 Trios

# Flag to process to L0 and L1 levels
PROC_L1=False

# Flag to process L1 to L2 levels
PROC_L2=True

# code path parameters:
if [ "$OSX" = "True" ] ; then
   octave='/opt/local/bin/octave-cli'
   DATA_PATH='/Volumes/Rivendell/AMT/Ben_data'
   CODE_PATH=$DATA_PATH'/code/'
   LOG_PATH=$DATA_PATH'/logs/'
   echo '==============================='
   echo  'Beware merge script issues... '
   echo  'HAVE YOU SET ulimit -n on OSX?'
   echo '==============================='
   sleep 10
else
   octave='/usr/local/bin/octave'
   DATA_PATH='/data/lazarev1/backup/cruise_data/AMT24/DallOlmo/HSAS/Satcon_extracted/Physical_units/' # I am not sure what this is for exactly
   CODE_PATH='/data/datasets/cruise_data/active/AMT24/DallOlmo/source/Radiometry/hsas.source/'
   LOG_PATH='/data/datasets/cruise_data/active/AMT24/DallOlmo/source/Radiometry/hsas.source/logs/'
fi

# Hsas parameters:
HSAS_CALIBRATED_PATH=$DATA_PATH # this is where the input calibrated HSAS data in SatCon format are 
HSAS_PROCESSED_PATH=$DATA_PATH'../../Processed_final' # this is where the processed data will go

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# OFF WE GO: no user configuration parameters below this point
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#make sure we run the last day
if [ "$OSX" = "True" ] ; then
   LAST_DATE=$(gdate -I -d "$LAST_DATE + 1 day")
else
   LAST_DATE=$(date -I -d "$LAST_DATE + 1 day")
fi

# L0....code removed as not needed for Hypersas



# L1....
if [ "$PROC_L1" = "True" ] ; then
  echo "L1 processing..."
  SCROLL_DATE=$FIRST_DATE

  while [ "$SCROLL_DATE" != $LAST_DATE ]; do
    # Loop through days
    if [ "$RUN_HSAS" = "True" ] ; then
      if [ "$OSX" = "True" ] ; then
        YYYY_DOY=$(gdate -d "$SCROLL_DATE" +%Y'-'%j)
        DOY=$(gdate -d "$SCROLL_DATE" +%j)
      else
        YYYY_DOY=$(date -d "$SCROLL_DATE" +%Y'-'%j)
        DOY=$(date -d "$SCROLL_DATE" +%j)
      fi
      if [ "$VERSION" = "v1" ] ; then
         echo 'Running with standard calibration'
         ls $HSAS_CALIBRATED_PATH/$DOY/*$YYYY_$DOY*H*.dat | head -n 1 | xargs -n1 $octave -qf proc_L0_L1.m hsas $DATA_PATH $VERSION #> $LOG_PATH/'Hsas_L1_'$SCROLL_DATE.'log'
      fi
    fi

    if [ "$OSX" = "True" ] ; then
       SCROLL_DATE=$(gdate -I -d "$SCROLL_DATE + 1 day")
    else
       SCROLL_DATE=$(date -I -d "$SCROLL_DATE + 1 day")
       echo $FIRST_DATE" "$SCROLL_DATE" "$LAST_DATE
    fi
  done
fi




# L2...
if [ "$PROC_L2" = "True" ] ; then
  echo "L2 processing..."
  SCROLL_DATE=$FIRST_DATE
  myDATE="${SCROLL_DATE//-/_}"

  while [ "$SCROLL_DATE" != $LAST_DATE ]; do
  # Loop through days
    if [ "$OSX" = "True" ] ; then
       DOY=$(gdate -d "$SCROLL_DATE" +%j)
    else
       DOY=$(date -d "$SCROLL_DATE" +%j)
    fi
    echo "---------------------------------------------------------------------"
    echo "L2 processing: "$DOY"/"$myDATE"_<INSTRUMENT>_L1_"$VERSION".mat"

    if [ "$RUN_TRIOS" = "True" ] ; then
      ls $TRIOS_PROCESSED_PATH/L1/$DOY/$myDATE'_Trios_L1_'$VERSION'.mat' | xargs -n1 octave -qf proc_L2_recal.m trios $DATA_PATH $VERSION 'cont' #> $LOG_PATH/'Trios_L2_'$SCROLL_DATE.'log' # both uses stricter filter for variability...needs work
    fi

    if [ "$RUN_HSAS" = "True" ] ; then

      # deprecated to * in new file nameing convention
      if [ "$OSX" = "True" ] ; then
        YYYY_DOY=$(gdate -d "$SCROLL_DATE" +%Y'-'%j)
        DOY=$(gdate -d "$SCROLL_DATE" +%j)
      else
        YYYY_DOY=$(date -d "$SCROLL_DATE" +%Y'-'%j)
        DOY=$(date -d "$SCROLL_DATE" +%j)
      fi
      # *

      if [ "$VERSION" = "v1" ] ; then
         echo 'HSAS: Running with standard calibration'
         ls $HSAS_PROCESSED_PATH/L1/$DOY/*$my_DATE'_L1_'$VERSION'.mat' | xargs -n1 -t octave -qf proc_L2.m hsas $DATA_PATH $VERSION 'continuous' #use continuous 
      fi
    fi

    if [ "$OSX" = "True" ] ; then
       SCROLL_DATE=$(gdate -I -d "$SCROLL_DATE + 1 day")
    else
       SCROLL_DATE=$(date -I -d "$SCROLL_DATE + 1 day")
    fi
    myDATE="${SCROLL_DATE//-/_}"
  done
fi
echo 'Done'
#-------------------------------------------------------------------------------

