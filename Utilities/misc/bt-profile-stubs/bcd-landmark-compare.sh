#!/usr/bin/env bash
# bt-profile-stubs/bcd-landmark-compare.sh
# Profile BCDTestForLandmarkCompare under macOS `sample`.
# This is the slowest BCD test (~1162 s) — dominant hotspot is the MSP
# exhaustive search in itkReflectiveCorrelationCenterToImageMetric.
#
# Usage (standalone):
#   bcd-landmark-compare.sh --build-dir DIR --data-dir DIR [--output-dir DIR]
#
# Usage (sourced by bt-profile-tests.sh):
#   source bcd-landmark-compare.sh   # exports _stub_name, _stub_cmd[], _stub_out
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

_stub_name="BCDTestForLandmarkCompare"
_stub_out="${OUTPUT_DIR}/${_stub_name}.sample.txt"

PROF_BCD_TS="${INNER_BUILD_DIR}/BRAINSConstellationDetector/TestSuite"
mkdir -p "${PROF_BCD_TS}"

_stub_cmd=(
  "${INNER_BUILD_DIR}/BRAINSConstellationDetector/TestSuite/BRAINSConstellationDetectorTestDriver"
  "BRAINSConstellationDetectorTest"
  "--inputVolume"    "${DATA_DIR}/T1.nii.gz"
  "--LLSModel"
    "${DATA_DIR}/Transforms_h5/LLSModel_50Lmks.h5"
  "--inputTemplateModel" "${DATA_DIR}/T1_50Lmks.mdl"
  "--outputVolume"
    "${PROF_BCD_TS}/BCDTestForLandmarkCompare_aligned.nrrd"
  "--outputLandmarksInACPCAlignedSpace"
    "${PROF_BCD_TS}/BCDTestForLandmarkCompare_ACPCSpace_aligned.fcsv"
  "--outputLandmarksInInputSpace"
    "${PROF_BCD_TS}/BCDTestForLandmarkCompare_InputSpace_aligned.fcsv"
  "--outputVerificationScript"
    "${PROF_BCD_TS}/BCDTestForLandmarkCompare_aligned.sh"
)

if [[ "${STANDALONE}" -eq 1 ]]; then
  "${SCRIPT_DIR}/../bt-profile-tests.sh" \
    --build-dir "${INNER_BUILD_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --stub "${_stub_name}"
fi
