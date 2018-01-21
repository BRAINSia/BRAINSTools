#!/bin/bash -
#===============================================================================
#
#          FILE: deface_registration_prep_t1.sh
#
#         USAGE: ./deface_registration_prep_t1.sh
#
#   DESCRIPTION:
## find /Shared/paulsen/Experiments/20160520_PREDICTHD_long_Results/*/*/*/TissueClassify/ -name t1_average_BRAINSABC.nii.gz > all_t1s.list
## for tt in $(cat all_t1s.list) ;do  bash ./deface_registration_prep_t1.sh $tt ; done
#
#         NOTES: This is a hard-coded reference script for aligning raw ACPC_RAW images to BAW outputs
#        AUTHOR: Hans J. Johnson
#  ORGANIZATION:
#       CREATED: 08/20/2017 07:24
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

export NSLOTS=2
export PATH=/Shared/pinc/sharedopt/20170302/RHEL7/NEP-intel/bin:$PATH
export LD_LIBRARY_PATH=/Shared/pinc/sharedopt/20170302/RHEL7/NEP-intel/bin:${LD_LIBRARY_PATH}

if [[ $# -ne 1 ]]; then
   echo "USING TEST CASE: $# arguments provided "
   IMG_ACPC=/Shared/paulsen/Experiments/20160520_PREDICTHD_long_Results/*/*/*/TissueClassify/t1_average_BRAINSABC.nii.gz
else
   IMG_ACPC=$1
fi

DEFACE_DIR=$(dirname $(dirname $IMG_ACPC))/DEFACE
mkdir -p ${DEFACE_DIR}

SITE=$(echo ${IMG_ACPC} |awk -F/ '{print $6}')
SUBJ=$(echo ${IMG_ACPC} |awk -F/ '{print $7}')
SESS=$(echo ${IMG_ACPC} |awk -F/ '{print $8}')

ORIG_DIR=/Shared/paulsen/MRx/${SITE}/${SUBJ}/${SESS}/ANONRAW

for img_orig in ${ORIG_DIR}/*T1-[13][50]_*.nii.gz; do
    ls -d ${img_orig}
    OUTDIR=${DEFACE_DIR}
    OUTFILE=$(basename ${img_orig//.nii.gz/_acpc.nii.gz})
    OUTVOL=${OUTDIR}/${OUTFILE}

if [[ ! -f ${OUTVOL} ]]; then
  echo "Creating: ${OUTVOL} "
  BRAINSFit \
    --costMetric "MSE" \
    --failureExitCode "-1" \
    --outputVolumePixelType "short" \
    --interpolationMode "ResampleInPlace" \
    --initializeTransformMode useCenterOfHeadAlign \
    --maskInferiorCutOffFromCenter  50 \
    --maskProcessingMode ROIAUTO  \
    --ROIAutoClosingSize 12 \
    --fixedVolume "${IMG_ACPC}" \
    --movingVolume "${img_orig}" \
    --outputVolume "${OUTVOL}" \
    --outputTransform "${OUTDIR}/BRAINSFitTest_RigidRotGeomNoMasks.mat" \
    --numberOfIterations 7500,7500,7500 \
    --numberOfHistogramBins 200  \
    --maximumStepLength 0.05  \
    --minimumStepLength "0.005,0.0005,0.000005" \
    --transformType "Rigid,Rigid,Rigid" \
    --relaxationFactor 0.5  \
    --translationScale 500  \
    --reproportionScale 1  \
     --histogramMatch \
    --skewScale 1
else
   echo "${OUTVOL} already created"
fi

done

#    --debugLevel 10
#    --numberOfIterations "7500" \
#    --numberOfHistogramBins "200" \
#    --translationScale "500" \
