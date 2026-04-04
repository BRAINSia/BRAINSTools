#!/usr/bin/env bash
# bt-profile-tests.sh
# Profile the slowest BRAINSFit and BRAINSABC tests using macOS `sample`.
#
# Usage:
#   bt-profile-tests.sh [OPTIONS]
#
# Options:
#   -b, --build-dir DIR   Inner build directory to use for test binaries
#                         (default: auto-detect Profiling build, fall back to Release)
#   -o, --output-dir DIR  Directory for sample output files  (default: /tmp/bt-profiles)
#   -h, --help            Show this message
#
# Requirements:
#   macOS only — uses the `sample` CLI tool bundled with Xcode Command Line Tools.
#   Build must have been made with BUILD_PROFILING=ON (or Debug) for accurate
#   symbol resolution.  The profiling inner build is at:
#     <superbuild>/BRAINSTools-Release-Profiling-EPRelease-build/
#
# Output:
#   One <testname>.sample.txt per test in OUTPUT_DIR.
#   A summary hotspot report printed to stdout at the end.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SUPERBUILD_DIR="${BRAINSTOOLS_BUILD_DIR:-${SOURCE_DIR}/build}"
OUTPUT_DIR="/tmp/bt-profiles"

# -----------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------
INNER_BUILD_DIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--build-dir) INNER_BUILD_DIR="$2"; shift 2 ;;
    -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,/^set -/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# Auto-detect inner build: prefer Profiling, fall back to Release
if [[ -z "${INNER_BUILD_DIR}" ]]; then
  for candidate in \
    "${SUPERBUILD_DIR}/BRAINSTools-Release-Profiling-EPRelease-build" \
    "${SUPERBUILD_DIR}/BRAINSTools-Release-EPRelease-build"
  do
    if [[ -d "${candidate}/bin" ]]; then
      INNER_BUILD_DIR="${candidate}"
      break
    fi
  done
fi

if [[ -z "${INNER_BUILD_DIR}" || ! -d "${INNER_BUILD_DIR}/bin" ]]; then
  echo "ERROR: Could not find an inner build directory. Use --build-dir or build first." >&2
  exit 1
fi

echo "=== BRAINSTools Performance Profiler ==="
echo "  Inner build : ${INNER_BUILD_DIR}"
echo "  Output dir  : ${OUTPUT_DIR}"
echo ""

mkdir -p "${OUTPUT_DIR}"

# Shared data directory (same for all inner builds)
DATA_DIR="${SUPERBUILD_DIR}/BRAINSTools-Release-EPRelease-build/ExternalData/TestData"
if [[ ! -d "${DATA_DIR}" ]]; then
  echo "ERROR: ExternalData not found at ${DATA_DIR}" >&2
  exit 1
fi

# Map paths from the Release build to the profiling build for output files
PROF_BRAINSFIT_TS="${INNER_BUILD_DIR}/BRAINSFit/TestSuite"
PROF_BRAINSABC_TS="${INNER_BUILD_DIR}/BRAINSABC/TestSuite"
mkdir -p "${PROF_BRAINSFIT_TS}" "${PROF_BRAINSABC_TS}"

# -----------------------------------------------------------------------
# Helper: run binary, sample it, report top hotspots
# -----------------------------------------------------------------------
_profile_test() {
  local test_name="$1"
  local sample_file="${OUTPUT_DIR}/${test_name}.sample.txt"
  shift
  local -a cmd=("$@")

  echo "--- Profiling: ${test_name} ---"
  echo "  Command: ${cmd[*]}"
  echo ""

  # Launch process in background
  "${cmd[@]}" &
  local pid=$!

  # Sample for up to 600 seconds (sample exits when PID ends)
  sample "${pid}" 600 -wait -f "${sample_file}" 2>/dev/null || true

  # Wait for the actual process to finish
  local status=0
  wait "${pid}" || status=$?

  if [[ ${status} -ne 0 ]]; then
    echo "  WARNING: ${test_name} exited with code ${status}"
  else
    echo "  PASSED: ${test_name}"
  fi

  echo "  Sample file: ${sample_file}"
  echo ""
}

# -----------------------------------------------------------------------
# Helper: print top N hotspot lines from a sample file
# -----------------------------------------------------------------------
_top_hotspots() {
  local label="$1"
  local sample_file="$2"
  local top_n="${3:-20}"

  echo "=== Top ${top_n} hotspots: ${label} ==="
  if [[ ! -f "${sample_file}" ]]; then
    echo "  (sample file not found: ${sample_file})"
    echo ""
    return
  fi

  # sample output format: leading spaces + count + spaces + symbol
  # Extract lines with numeric counts, sort by count descending, deduplicate symbols
  grep -E '^ +[0-9]+ +[A-Za-z_]' "${sample_file}" \
    | awk '{count=$1; $1=""; sym=$0; gsub(/^ +/, "", sym); print count, sym}' \
    | sort -rn \
    | awk '!seen[$2]++' \
    | head -"${top_n}" \
    | awk '{printf "  %6d  %s\n", $1, $2}' \
    || echo "  (no symbols found — check sample file)"

  echo ""
}

# -----------------------------------------------------------------------
# Test 1: BRAINSABCSmallTest  (dominant, ~140s)
# -----------------------------------------------------------------------
_profile_test "BRAINSABCSmallTest" \
  "${INNER_BUILD_DIR}/BRAINSABC/TestSuite/BRAINSABCTestDriver" \
  "--compare" \
    "${DATA_DIR}/BRAINSABCSmallLabels.nii.gz" \
    "${PROF_BRAINSABC_TS}/BRAINSABCSmallLabels.test.nii.gz" \
  "--compareIntensityTolerance" "1" \
  "--compareRadiusTolerance" "1" \
  "--compareNumberOfPixelsTolerance" "10000" \
  "BRAINSABCTest" \
  "--atlasDefinition" \
    "${INNER_BUILD_DIR}/BRAINSABC/TestSuite/BRAINSABCSmallExtendedAtlasDefinition.xml" \
  "--atlasToSubjectInitialTransform" \
    "${DATA_DIR}/BRAINSABCSmall_atlas_to_subject_transform.h5" \
  "--atlasToSubjectTransform" "BRAINSABCSmall_atlas_to_subject_transform.h5" \
  "--atlasToSubjectTransformType" "Affine" \
  "--debuglevel" "0" \
  "--filterIteration" "0" \
  "--filterMethod" "GradientAnisotropicDiffusion" \
  "--gridSize" "10,10,10" \
  "--inputVolumeTypes" "T1,T2" \
  "--inputVolumes" "${DATA_DIR}/affine_t1.nrrd" \
  "--inputVolumes" "${DATA_DIR}/affine_t2.nrrd" \
  "--interpolationMode" "Linear" \
  "--maxBiasDegree" "4" \
  "--maxIterations" "1" \
  "--outputDir" "${PROF_BRAINSABC_TS}/" \
  "--outputDirtyLabels" "${PROF_BRAINSABC_TS}/BRAINSABCSmallvolume_label_seg.nii.gz" \
  "--outputFormat" "NIFTI" \
  "--outputLabels" "${PROF_BRAINSABC_TS}/BRAINSABCSmallLabels.test.nii.gz" \
  "--outputVolumes" "${PROF_BRAINSABC_TS}/BRAINSABCSmallT1_1.nii.gz" \
  "--outputVolumes" "${PROF_BRAINSABC_TS}/BRAINSABCSmallT2_1.nii.gz" \
  "--posteriorTemplate" "${PROF_BRAINSABC_TS}/BRAINSABCSmallPOST_%s.nii.gz" \
  "--purePlugsThreshold" "0.2"

# -----------------------------------------------------------------------
# Test 2: BSplineOnlyRescaleHeadMasks  (~16s)
# -----------------------------------------------------------------------
_profile_test "BRAINSFitTest_BSplineOnlyRescaleHeadMasks" \
  "${INNER_BUILD_DIR}/BRAINSFit/TestSuite/BRAINSFitTestDriver" \
  "--compare" \
    "${DATA_DIR}/BRAINSFitTest_BSplineOnlyRescaleHeadMasks.result.nii.gz" \
    "${PROF_BRAINSFIT_TS}/BRAINSFitTest_BSplineOnlyRescaleHeadMasks.test.nii.gz" \
  "--compareIntensityTolerance" "9" \
  "--compareRadiusTolerance" "1" \
  "--compareNumberOfPixelsTolerance" "300" \
  "BRAINSFitTest" \
  "--costMetric" "MMI" \
  "--failureExitCode" "-1" \
  "--writeTransformOnFailure" \
  "--numberOfIterations" "1500" \
  "--numberOfHistogramBins" "200" \
  "--splineGridSize" "7,5,6" \
  "--samplingPercentage" "0.5" \
  "--translationScale" "250" \
  "--minimumStepLength" "0.01" \
  "--outputVolumePixelType" "short" \
  "--maskProcessingMode" "ROIAUTO" \
  "--initialTransform" \
    "${DATA_DIR}/Transforms_h5/Initializer_BRAINSFitTest_BSplineAnteScaleRotationRescaleHeadMasks.h5" \
  "--transformType" "BSpline" \
  "--fixedVolume" "${DATA_DIR}/test.nii.gz" \
  "--movingVolume" "${DATA_DIR}/rotation.rescale.rigid.nii.gz" \
  "--outputVolume" "${PROF_BRAINSFIT_TS}/BRAINSFitTest_BSplineOnlyRescaleHeadMasks.test.nii.gz" \
  "--outputTransform" "${PROF_BRAINSFIT_TS}/BRAINSFitTest_BSplineOnlyRescaleHeadMasks.h5" \
  "--debugLevel" "10" \
  "--maxBSplineDisplacement" "7.3" \
  "--projectedGradientTolerance" "1e-4" \
  "--costFunctionConvergenceFactor" "1e+9"

# -----------------------------------------------------------------------
# Test 3: MIHAffineRotationMasks  (~12s)
# -----------------------------------------------------------------------
_profile_test "BRAINSFitTest_MIHAffineRotationMasks" \
  "${INNER_BUILD_DIR}/BRAINSFit/TestSuite/BRAINSFitTestDriver" \
  "--compare" \
    "${DATA_DIR}/BRAINSFitTest_MIHAffineRotationMasks.result.nii.gz" \
    "${PROF_BRAINSFIT_TS}/BRAINSFitTest_MIHAffineRotationMasks.test.nii.gz" \
  "--compareIntensityTolerance" "7" \
  "--compareRadiusTolerance" "0" \
  "--compareNumberOfPixelsTolerance" "777" \
  "BRAINSFitTest" \
  "--costMetric" "MIH" \
  "--failureExitCode" "-1" \
  "--writeTransformOnFailure" \
  "--numberOfIterations" "2500" \
  "--numberOfHistogramBins" "200" \
  "--samplingPercentage" "0.5" \
  "--translationScale" "250" \
  "--minimumStepLength" "0.001" \
  "--outputVolumePixelType" "uchar" \
  "--transformType" "Affine" \
  "--initialTransform" \
    "${DATA_DIR}/Transforms_h5/BRAINSFitTest_Initializer_RigidRotationNoMasks.h5" \
  "--maskProcessingMode" "ROI" \
  "--fixedVolume" "${DATA_DIR}/test.nii.gz" \
  "--fixedBinaryVolume" "${DATA_DIR}/test_mask.nii.gz" \
  "--movingVolume" "${DATA_DIR}/rotation.test.nii.gz" \
  "--movingBinaryVolume" "${DATA_DIR}/rotation.test_mask.nii.gz" \
  "--outputVolume" "${PROF_BRAINSFIT_TS}/BRAINSFitTest_MIHAffineRotationMasks.test.nii.gz" \
  "--outputTransform" "${PROF_BRAINSFIT_TS}/BRAINSFitTest_MIHAffineRotationMasks.h5" \
  "--debugLevel" "50"

# -----------------------------------------------------------------------
# Hotspot Report
# -----------------------------------------------------------------------
echo ""
echo "================================================================"
echo " HOTSPOT REPORT"
echo "================================================================"
echo ""

_top_hotspots "BRAINSABCSmallTest"                    "${OUTPUT_DIR}/BRAINSABCSmallTest.sample.txt"
_top_hotspots "BRAINSFit BSplineOnlyRescaleHeadMasks"  "${OUTPUT_DIR}/BRAINSFitTest_BSplineOnlyRescaleHeadMasks.sample.txt"
_top_hotspots "BRAINSFit MIHAffineRotationMasks"       "${OUTPUT_DIR}/BRAINSFitTest_MIHAffineRotationMasks.sample.txt"

echo "Sample files saved to: ${OUTPUT_DIR}/"
echo ""
echo "To explore interactively, open in Instruments:"
echo "  xcrun xctrace import --input <file> --output out.trace && open out.trace"
