#!/usr/bin/env bash
# bt-profile-stubs/brainsabc-small.sh
# Profile BRAINSABCSmallTest under macOS `sample`.
#
# Usage (standalone):
#   brainsabc-small.sh --build-dir DIR --data-dir DIR [--output-dir DIR]
#
# Usage (sourced by bt-profile-tests.sh):
#   source brainsabc-small.sh   # exports _stub_name, _stub_cmd[], _stub_out
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

_stub_name="BRAINSABCSmallTest"
_stub_out="${OUTPUT_DIR}/${_stub_name}.sample.txt"

PROF_BRAINSABC_TS="${INNER_BUILD_DIR}/BRAINSABC/TestSuite"
mkdir -p "${PROF_BRAINSABC_TS}"

_stub_cmd=(
  "${INNER_BUILD_DIR}/BRAINSABC/TestSuite/BRAINSABCTestDriver"
  "--compare"
    "${DATA_DIR}/BRAINSABCSmallLabels.nii.gz"
    "${PROF_BRAINSABC_TS}/BRAINSABCSmallLabels.test.nii.gz"
  "--compareIntensityTolerance" "1"
  "--compareRadiusTolerance"    "1"
  "--compareNumberOfPixelsTolerance" "10000"
  "BRAINSABCTest"
  "--atlasDefinition"
    "${INNER_BUILD_DIR}/BRAINSABC/TestSuite/BRAINSABCSmallExtendedAtlasDefinition.xml"
  "--atlasToSubjectInitialTransform"
    "${DATA_DIR}/BRAINSABCSmall_atlas_to_subject_transform.h5"
  "--atlasToSubjectTransform" "BRAINSABCSmall_atlas_to_subject_transform.h5"
  "--atlasToSubjectTransformType" "Affine"
  "--debuglevel" "0"
  "--filterIteration" "0"
  "--filterMethod" "GradientAnisotropicDiffusion"
  "--gridSize" "10,10,10"
  "--inputVolumeTypes" "T1,T2"
  "--inputVolumes" "${DATA_DIR}/affine_t1.nrrd"
  "--inputVolumes" "${DATA_DIR}/affine_t2.nrrd"
  "--interpolationMode" "Linear"
  "--maxBiasDegree" "4"
  "--maxIterations" "1"
  "--outputDir" "${PROF_BRAINSABC_TS}/"
  "--outputDirtyLabels"
    "${PROF_BRAINSABC_TS}/BRAINSABCSmallvolume_label_seg.nii.gz"
  "--outputFormat" "NIFTI"
  "--outputLabels"
    "${PROF_BRAINSABC_TS}/BRAINSABCSmallLabels.test.nii.gz"
  "--outputVolumes" "${PROF_BRAINSABC_TS}/BRAINSABCSmallT1_1.nii.gz"
  "--outputVolumes" "${PROF_BRAINSABC_TS}/BRAINSABCSmallT2_1.nii.gz"
  "--posteriorTemplate"
    "${PROF_BRAINSABC_TS}/BRAINSABCSmallPOST_%s.nii.gz"
  "--purePlugsThreshold" "0.2"
)

if [[ "${STANDALONE}" -eq 1 ]]; then
  # Run the stub directly, then print hotspots
  "${SCRIPT_DIR}/../bt-profile-tests.sh" \
    --build-dir "${INNER_BUILD_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --stub "${_stub_name}"
fi
