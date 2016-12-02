#!/bin/bash -ex

SLICER_EXTENSION_BIN=/Users/johnsonhj/Desktop/Slicer.app/Contents/Extensions-25603/SlicerDMRI/lib/Slicer-4.7/cli-modules

DWIIN=$1
BASENAME=${DWIIN//.nhdr/}
DTIOUT=${BASENAME}_SlicerDTI.nrrd
DWIBRAINMASK=${BASENAME}_BrainMask.nrrd
BASELINEOUT=${BASENAME}_Baseline.nrrd
BASELINEOUTDTI=${BASENAME}_DTIBaseline.nrrd
FA=${BASENAME}_SlicerFA.nrrd
MD=${BASENAME}_SlicerMD.nrrd



#<std::string> (required)  Input DWI volume
#<std::string> (required)  Extracted baseline volume
#<std::string> (required)  Output Diffusion Brain Mask
${SLICER_EXTENSION_BIN}/DiffusionWeightedVolumeMasking \
  --baselineBValueThreshold 100 \
  --removeislands \
  ${DWIIN} ${BASELINEOUT} ${DWIBRAINMASK}

#   <std::string> (required)  Input Diffusion Weighted Image (DWI) volume
#   <std::string> (required)  Estimated Diffusion Tensor Image (DTI) volume
#   <std::string> (required)  Estimated baseline (non-Diffusion Weighted) volume
${SLICER_EXTENSION_BIN}/DWIToDTIEstimation \
  --enumeration LS \
  --mask ${DWIBRAINMASK} \
  ${DWIIN} ${DTIOUT} ${BASELINEOUTDTI}


#  -e <FractionalAnisotropy|Trace|Determinant|RelativeAnisotropy|Mode
#     |LinearMeasure|PlanarMeasure|SphericalMeasure|MinEigenvalue
#     |MidEigenvalue|MaxEigenvalue|ParallelDiffusivity
#     |PerpendicularDiffusivity>,  --enumeration <FractionalAnisotropy
#     |Trace|Determinant|RelativeAnisotropy|Mode|LinearMeasure
#     |PlanarMeasure|SphericalMeasure|MinEigenvalue|MidEigenvalue
#     |MaxEigenvalue|ParallelDiffusivity|PerpendicularDiffusivity>
#    Type of scalar measurement to perform (default: FractionalAnisotropy)

# <std::string> (required)  Input DTI volume
# <std::string> (required)  Scalar volume derived from tensor
${SLICER_EXTENSION_BIN}/DiffusionTensorScalarMeasurements \
  --enumeration FractionalAnisotropy \
  ${DTIOUT} \
  ${FA}


##=======================
if [ 1 -eq 1 ]; then
  #  --estimationMethod <lls|wls|nls|ml>
  #DTIESTIMBIN="/Users/johnsonhj/src/NEP-11/bin/dtiestim"
  DTIESTIMBIN="/Users/johnsonhj/src/DTIProcessToolkit-build/DTIProcess-build/bin/dtiestim"
  DTIPROCESSBIN="/Users/johnsonhj/src/DTIProcessToolkit-build/DTIProcess-build/bin/dtiprocess"

  ${DTIESTIMBIN}  --estimationMethod wls  \
     --inputBrainMaskVolume ${DWIBRAINMASK} \
     --inputDWIVolume ${DWIIN} \
     --outputDTIVolume ${BASENAME}_ESTIMDTI.nhdr
  ${DTIPROCESSBIN} --inputDTIVolume ${BASENAME}_ESTIMDTI.nhdr \
        --outputFAVolume ${BASENAME}_PROCFA.nrrd \
        --outputMDVolume ${BASENAME}_PROCMD.nrrd \
        --saveScalarsAsFloat

fi
