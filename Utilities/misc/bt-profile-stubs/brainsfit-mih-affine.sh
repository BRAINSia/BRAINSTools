#!/usr/bin/env bash
# bt-profile-stubs/brainsfit-mih-affine.sh
# Profile BRAINSFitTest_MIHAffineRotationMasks under macOS `sample`.
#
# Usage (standalone):
#   brainsfit-mih-affine.sh --build-dir DIR --data-dir DIR [--output-dir DIR]
#
# Usage (sourced by bt-profile-tests.sh):
#   source brainsfit-mih-affine.sh   # exports _stub_name, _stub_cmd[], _stub_out
#
# Required variables when sourced:
#   INNER_BUILD_DIR  DATA_DIR  OUTPUT_DIR
set -euo pipefail

# ---------- resolve args when run standalone --------------------------------
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  INNER_BUILD_DIR="" DATA_DIR="" OUTPUT_DIR="/tmp/bt-profiles"
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --build-dir)  INNER_BUILD_DIR="$2"; shift 2 ;;
      --data-dir)   DATA_DIR="$2";        shift 2 ;;
      --output-dir) OUTPUT_DIR="$2";      shift 2 ;;
      *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
  done
  if [[ -z "${INNER_BUILD_DIR}" || -z "${DATA_DIR}" ]]; then
    echo "Usage: $0 --build-dir DIR --data-dir DIR [--output-dir DIR]" >&2
    exit 1
  fi
  STANDALONE=1
else
  STANDALONE=0
fi
# ---------------------------------------------------------------------------

_stub_name="BRAINSFitTest_MIHAffineRotationMasks"
_stub_out="${OUTPUT_DIR}/${_stub_name}.sample.txt"

PROF_BRAINSFIT_TS="${INNER_BUILD_DIR}/BRAINSFit/TestSuite"
mkdir -p "${PROF_BRAINSFIT_TS}"

_stub_cmd=(
  "${INNER_BUILD_DIR}/BRAINSFit/TestSuite/BRAINSFitTestDriver"
  "--compare"
    "${DATA_DIR}/BRAINSFitTest_MIHAffineRotationMasks.result.nii.gz"
    "${PROF_BRAINSFIT_TS}/BRAINSFitTest_MIHAffineRotationMasks.test.nii.gz"
  "--compareIntensityTolerance" "7"
  "--compareRadiusTolerance"    "0"
  "--compareNumberOfPixelsTolerance" "777"
  "BRAINSFitTest"
  "--costMetric"       "MIH"
  "--failureExitCode"  "-1"
  "--writeTransformOnFailure"
  "--numberOfIterations"    "2500"
  "--numberOfHistogramBins" "200"
  "--samplingPercentage"    "0.5"
  "--translationScale"      "250"
  "--minimumStepLength"     "0.001"
  "--outputVolumePixelType" "uchar"
  "--transformType"  "Affine"
  "--initialTransform"
    "${DATA_DIR}/Transforms_h5/BRAINSFitTest_Initializer_RigidRotationNoMasks.h5"
  "--maskProcessingMode"   "ROI"
  "--fixedVolume"          "${DATA_DIR}/test.nii.gz"
  "--fixedBinaryVolume"    "${DATA_DIR}/test_mask.nii.gz"
  "--movingVolume"         "${DATA_DIR}/rotation.test.nii.gz"
  "--movingBinaryVolume"   "${DATA_DIR}/rotation.test_mask.nii.gz"
  "--outputVolume"
    "${PROF_BRAINSFIT_TS}/BRAINSFitTest_MIHAffineRotationMasks.test.nii.gz"
  "--outputTransform"
    "${PROF_BRAINSFIT_TS}/BRAINSFitTest_MIHAffineRotationMasks.h5"
  "--debugLevel" "0"
)

if [[ "${STANDALONE}" -eq 1 ]]; then
  "${SCRIPT_DIR}/../bt-profile-tests.sh" \
    --build-dir "${INNER_BUILD_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --stub "${_stub_name}"
fi
