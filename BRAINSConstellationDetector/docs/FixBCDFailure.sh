#/bin/bash
#
# \author Hans J. Johnson
#
# This script is designed to assist with fixing failed BCD
# data processing.


while getopts "c:p:" opt; do
  case $opt in
     p)
	echo "Processing phase 1"
	PHASE_TO_RUN=$OPTARG
        ;;
     c)
        CONFIG_FILE=$OPTARG
        echo "For configuration ${CONFIG_FILE}"
        ;;
     *)
	echo "ERROR:  Invalid option"
        exit -1;
        ;;
   esac
done


if [[ ! -f ${CONFIG_FILE} ]] ;then
   echo "ERROR: Can not source ${CONFIG_FILE}"
   exit -1;
fi

source ${CONFIG_FILE}

MANUAL_NOT_DONE_SENTINAL=MANUAL_LANDMARKS_NOT_SET.txt

mkdir -p ${OUTPUT_DIR}
cp $0 ${OUTPUT_DIR}
cp ${CONFIG_FILE} ${OUTPUT_DIR}
pushd ${OUTPUT_DIR}



##
#
# Phase 1  (Implicit output of EMSP.nrrd)
PHASE1_OUTPUT_VOL=${OUTPUT_FILENAME}
if [[ ${PHASE_TO_RUN} -eq 1 ]]; then
  if [[ -f EMSP.nrrd ]] || [[ -f EMSP.fcsv ]] ;then
     echo "ERROR:  files exist already.  Remove EMSP.nrrd and EMSP.fcsv before continuing."
     exit -1;
  fi
  ${PROG_PATH}/BRAINSConstellationDetector \
    --inputVolume ${INPUT_VOL} \
    --acLowerBound 80.000000 \
    --LLSModel ${LLS_MODEL} \
    --inputTemplateModel ${TEMP_MODEL} \
    --houghEyeDetectorMode 1 \
    --interpolationMode Linear \
    --forceHoughEyeDetectorReportFailure

  echo "DELETE ME AFTER MANUAL EDITING OF LANDMARKS" > ${MANUAL_NOT_DONE_SENTINAL}
  echo "\n\n\n\n  Use slicer to edit the EMSP.nrrd and EMSP.fcsv files to manually"
  echo "  fix the eye locations\n\n"

  echo "Start Slicer3D"
  echo "Get data from $(pwd)"
  echo "Then delete the file ${MANUAL_NOT_DONE_SENTINAL}"

  mv EMSP.nrrd ${PHASE1_OUTPUT_VOL}
  touch ${PHASE1_OUTPUT_VOL}_noDenoise
fi


##
#
# Phase 2
PHASE_2_OUTPUT_VOL=step2_output.nii.gz
if [[ ${PHASE_TO_RUN} -eq 2 ]]; then
  ## NOTE THIS IS OFTEN NOT NECESSARY
  if [[ -f ${MANUAL_NOT_DONE_SENTINAL} ]]; then
     echo "ERROR: Complete manual editing of eye landmarks"
     echo "       then remove ${MANUAL_NOT_DONE_SENTINAL}"
     exit -1
  fi
  ${PROG_PATH}/BRAINSConstellationDetector \
    --inputVolume ${PHASE1_OUTPUT_VOL} \
    --acLowerBound 80.000000 \
    --LLSModel ${LLS_MODEL} \
    --inputTemplateModel ${TEMP_MODEL} \
    --houghEyeDetectorMode 1 \
    --interpolationMode Linear \
    --atlasLandmarkWeights ${ATLAS_LMKS_WTS} \
    --atlasLandmarks ${ATLAS_LMKS} \
    --atlasVolume ${ATLAS_VOLUME} \
    --outputResampledVolume ${PHASE_2_OUTPUT_VOL}
fi

##
#
# Phase 3
if [[ ${PHASE_TO_RUN} -eq 3 ]]; then
  ## NOTE THIS IS OFTEN NOT NECESSARY
  if [[ ! -f ${PHASE_2_OUTPUT_VOL} ]] ; then
    echo "ERROR:  You must complte phase 2 before processing phase 3"
    exit -1
  fi
  ${PROG_PATH}/BRAINSConstellationDetector \
    --inputVolume ${PHASE_2_OUTPUT_VOL} \
    --acLowerBound 80.000000 \
    --LLSModel ${LLS_MODEL} \
    --houghEyeDetectorMode 1 \
    --interpolationMode Linear \
    --inputTemplateModel ${TEMP_MODEL} \
    --atlasLandmarkWeights ${ATLAS_LMKS_WTS} \
    --atlasLandmarks ${ATLAS_LMKS} \
    --atlasVolume ${ATLAS_VOLUME} \
    --outputLandmarksInACPCAlignedSpace BCD_ACPC_Landmarks.fcsv \
    --outputLandmarksInInputSpace BCD_Original.fcsv \
    --outputResampledVolume BCD_ACPC.nii.gz

    echo "NOTE:  Remember to update input configuration files for BAW if necessary"
    echo "       /Shared/johnsonhj/TrackOn/reports/etc/latest_black_list_images.lst"
    echo "       edited*.csv"
    echo "       etc...."
fi
