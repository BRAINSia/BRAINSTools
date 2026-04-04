#!/usr/bin/env bash
# bt-profile-stubs/brainsfit-bspline-rescale.sh
# Profile BRAINSFitTest_BSplineOnlyRescaleHeadMasks under macOS `sample`.
#
# Usage (standalone):
#   brainsfit-bspline-rescale.sh --build-dir DIR --data-dir DIR [--output-dir DIR]
#
# Usage (sourced by bt-profile-tests.sh):
#   source brainsfit-bspline-rescale.sh   # exports _stub_name, _stub_cmd[], _stub_out
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

_stub_name="BRAINSFitTest_BSplineOnlyRescaleHeadMasks"
_stub_out="${OUTPUT_DIR}/${_stub_name}.sample.txt"

PROF_BRAINSFIT_TS="${INNER_BUILD_DIR}/BRAINSFit/TestSuite"
mkdir -p "${PROF_BRAINSFIT_TS}"

_stub_cmd=(
  "${INNER_BUILD_DIR}/BRAINSFit/TestSuite/BRAINSFitTestDriver"
  "--compare"
    "${DATA_DIR}/BRAINSFitTest_BSplineOnlyRescaleHeadMasks.result.nii.gz"
    "${PROF_BRAINSFIT_TS}/BRAINSFitTest_BSplineOnlyRescaleHeadMasks.test.nii.gz"
  "--compareIntensityTolerance" "9"
  "--compareRadiusTolerance"    "1"
  "--compareNumberOfPixelsTolerance" "300"
  "BRAINSFitTest"
  "--costMetric" "MMI"
  "--failureExitCode" "-1"
  "--writeTransformOnFailure"
  "--numberOfIterations"    "1500"
  "--numberOfHistogramBins" "100"
  "--splineGridSize"        "7,5,6"
  "--samplingPercentage"    "0.5"
  "--translationScale"      "250"
  "--minimumStepLength"     "0.01"
  "--outputVolumePixelType" "short"
  "--maskProcessingMode"    "ROIAUTO"
  "--initialTransform"
    "${DATA_DIR}/Transforms_h5/Initializer_BRAINSFitTest_BSplineAnteScaleRotationRescaleHeadMasks.h5"
  "--transformType" "BSpline"
  "--fixedVolume"   "${DATA_DIR}/test.nii.gz"
  "--movingVolume"  "${DATA_DIR}/rotation.rescale.rigid.nii.gz"
  "--outputVolume"
    "${PROF_BRAINSFIT_TS}/BRAINSFitTest_BSplineOnlyRescaleHeadMasks.test.nii.gz"
  "--outputTransform"
    "${PROF_BRAINSFIT_TS}/BRAINSFitTest_BSplineOnlyRescaleHeadMasks.h5"
  "--debugLevel"                    "10"
  "--maxBSplineDisplacement"        "7.3"
  "--projectedGradientTolerance"    "1e-4"
  "--costFunctionConvergenceFactor" "1e+9"
)

if [[ "${STANDALONE}" -eq 1 ]]; then
  "${SCRIPT_DIR}/../bt-profile-tests.sh" \
    --build-dir "${INNER_BUILD_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --stub "${_stub_name}"
fi
