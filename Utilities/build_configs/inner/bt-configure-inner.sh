#!/usr/bin/env bash
# bt-configure-inner.sh
# Configure a BRAINSTools inner (Phase-II) build, reusing ExternalProject
# trees that were built during the SuperBuild (Phase-I).
#
# The SuperBuild always builds external dependencies in Release mode
# (EXTERNAL_PROJECT_BUILD_TYPE=Release), so both a Debug and a Release
# inner build share the same EP trees.
#
# Usage:
#   bt-configure-inner.sh [OPTIONS]
#
# Options:
#   -t, --type TYPE        Build type: Debug|Release|RelWithDebInfo  (default: Debug)
#   -b, --build-dir DIR    SuperBuild binary directory               (default: auto-detect)
#   -s, --source-dir DIR   BRAINSTools source root                   (default: auto-detect)
#   -m, --modules FILE     CMake init-cache with USE_* flags         (default: sibling BRAINSTools-modules.cmake)
#   -f, --force            Re-configure even if CMakeCache.txt exists
#   -h, --help             Show this message
#
# Auto-detection:
#   SOURCE_DIR  — resolved from this script's location (../../.. from Utilities/build_configs/inner/)
#   BUILD_DIR   — $SOURCE_DIR/build if it contains a CMakeCache.txt (superbuild was run there)
#
# Environment:
#   BRAINSTOOLS_BUILD_DIR  overrides --build-dir auto-detection
#   BRAINSTOOLS_SOURCE_DIR overrides --source-dir auto-detection
#   PARALLEL_JOBS          passed as -j to ninja (default: nproc/2)

set -euo pipefail

# -----------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------
BUILD_TYPE="Debug"
FORCE=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source dir: three levels up from Utilities/build_configs/inner/
AUTO_SOURCE_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
SOURCE_DIR="${BRAINSTOOLS_SOURCE_DIR:-${AUTO_SOURCE_DIR}}"

# Build dir: $SOURCE_DIR/build (or override)
AUTO_BUILD_DIR="${SOURCE_DIR}/build"
SUPERBUILD_DIR="${BRAINSTOOLS_BUILD_DIR:-${AUTO_BUILD_DIR}}"

# Modules init-cache file
MODULES_FILE="${SCRIPT_DIR}/BRAINSTools-modules.cmake"

# -----------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--type)         BUILD_TYPE="$2"; shift 2 ;;
    -b|--build-dir)    SUPERBUILD_DIR="$2"; shift 2 ;;
    -s|--source-dir)   SOURCE_DIR="$2"; shift 2 ;;
    -m|--modules)      MODULES_FILE="$2"; shift 2 ;;
    -f|--force)        FORCE=1; shift ;;
    -h|--help)
      sed -n '2,/^set -/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

EP_BUILD_TYPE="Release"   # EPs are always built Release

# -----------------------------------------------------------------------
# Validate inputs
# -----------------------------------------------------------------------
if [[ ! -f "${SUPERBUILD_DIR}/CMakeCache.txt" ]]; then
  echo "ERROR: SuperBuild CMakeCache.txt not found at: ${SUPERBUILD_DIR}/CMakeCache.txt" >&2
  echo "       Run the SuperBuild first, or specify --build-dir." >&2
  exit 1
fi

if [[ ! -f "${MODULES_FILE}" ]]; then
  echo "ERROR: Modules init-cache not found: ${MODULES_FILE}" >&2
  exit 1
fi

# -----------------------------------------------------------------------
# Auto-detect EP paths from the SuperBuild CMakeCache.txt
# -----------------------------------------------------------------------
sb_cache="${SUPERBUILD_DIR}/CMakeCache.txt"

_cache_val() {
  grep "^${1}:" "${sb_cache}" 2>/dev/null | head -1 | cut -d= -f2
}

INSTALL_PREFIX="$(_cache_val CMAKE_INSTALL_PREFIX)"
if [[ -z "${INSTALL_PREFIX}" ]]; then
  echo "ERROR: Could not read CMAKE_INSTALL_PREFIX from ${sb_cache}" >&2
  exit 1
fi

# ITK cmake dir (glob for version-stamped directory)
ITK_CMAKE_DIR=$(find "${SUPERBUILD_DIR}/ITKv5-${EP_BUILD_TYPE}-build/lib/cmake" \
  -maxdepth 1 -name "ITK-*" -type d 2>/dev/null | sort | tail -1)
if [[ -z "${ITK_CMAKE_DIR}" ]]; then
  echo "ERROR: Could not find ITK cmake config under ${SUPERBUILD_DIR}/ITKv5-${EP_BUILD_TYPE}-build" >&2
  exit 1
fi

VTK_DIR="${SUPERBUILD_DIR}/VTK-${EP_BUILD_TYPE}-build"
SEM_DIR="${SUPERBUILD_DIR}/SlicerExecutionModel-${EP_BUILD_TYPE}-build"
TBB_DIR="${SUPERBUILD_DIR}/tbb-${EP_BUILD_TYPE}-build"
ANTS_SRC="${SUPERBUILD_DIR}/ANTs"
ANTS_LIB="${INSTALL_PREFIX}/lib"   # ANTs installed into the superbuild install prefix

for _dir in "${VTK_DIR}" "${SEM_DIR}" "${TBB_DIR}" "${ANTS_SRC}" "${ANTS_LIB}"; do
  if [[ ! -d "${_dir}" ]]; then
    echo "ERROR: Required EP directory not found: ${_dir}" >&2
    exit 1
  fi
done

# -----------------------------------------------------------------------
# Inner build directory
# -----------------------------------------------------------------------
INNER_BUILD_DIR="${SUPERBUILD_DIR}/BRAINSTools-${BUILD_TYPE}-EP${EP_BUILD_TYPE}-build"

if [[ -f "${INNER_BUILD_DIR}/CMakeCache.txt" && "${FORCE}" -eq 0 ]]; then
  echo "INFO: ${INNER_BUILD_DIR}/CMakeCache.txt already exists — skipping configure."
  echo "      Use --force to reconfigure."
  exit 0
fi

# -----------------------------------------------------------------------
# Optimization flags — match BRAINSTools-Refactor-Release.cmake for Release,
# add frame-pointer/sibling-call flags for Debug (profiler-friendly).
# These are consumed by Common.cmake which appends them to CMAKE_CXX_FLAGS.
# -----------------------------------------------------------------------
BASE_OPT_FLAGS="-mtune=native;-march=native;\
-Wtautological-overlap-compare;-Wtautological-compare;\
-Wtautological-bitwise-compare;-Wbitwise-conditional-parentheses;\
-Wrange-loop-analysis;-Wmisleading-indentation;\
-Wc99-designator;-Wreorder-init-list;\
-Wsizeof-array-div;-Wxor-used-as-pow;\
-Wfinal-dtor-non-final-class"

if [[ "${BUILD_TYPE}" == "Debug" ]]; then
  OPT_FLAGS="-fno-omit-frame-pointer;-fno-optimize-sibling-calls;${BASE_OPT_FLAGS}"
else
  OPT_FLAGS="${BASE_OPT_FLAGS}"
fi

# -----------------------------------------------------------------------
# Configure
# -----------------------------------------------------------------------
echo "=== Configuring BRAINSTools inner build ==="
echo "  BUILD_TYPE   : ${BUILD_TYPE}"
echo "  SOURCE_DIR   : ${SOURCE_DIR}"
echo "  INNER_BUILD  : ${INNER_BUILD_DIR}"
echo "  ITK_DIR      : ${ITK_CMAKE_DIR}"
echo "  VTK_DIR      : ${VTK_DIR}"
echo "  INSTALL_PFX  : ${INSTALL_PREFIX}"
echo ""

cmake \
  -C "${MODULES_FILE}" \
  -S "${SOURCE_DIR}" \
  -B "${INNER_BUILD_DIR}" \
  -G Ninja \
  -DCMAKE_BUILD_TYPE:STRING="${BUILD_TYPE}" \
  -DCMAKE_CXX_STANDARD:STRING=17 \
  -DLOCAL_PROJECT_NAME:STRING=BRAINSTools \
  -DSUPERBUILD_TOPLEVEL_PROJECT:STRING=BRAINSTools \
  -DBRAINSTools_SUPERBUILD:BOOL=OFF \
  -DITK_DIR:PATH="${ITK_CMAKE_DIR}" \
  -DVTK_DIR:PATH="${VTK_DIR}" \
  -DSlicerExecutionModel_DIR:PATH="${SEM_DIR}" \
  -DTBB_DIR:PATH="${TBB_DIR}" \
  -DZLIB_ROOT:PATH="${INSTALL_PREFIX}" \
  -DANTs_SOURCE_DIR:PATH="${ANTS_SRC}" \
  -DANTs_LIBRARY_DIR:PATH="${ANTS_LIB}" \
  -DBRAINSToools_CXX_OPTIMIZATION_FLAGS:STRING="${OPT_FLAGS}" \
  -DBRAINSToools_C_OPTIMIZATION_FLAGS:STRING="${OPT_FLAGS}"

echo ""
echo "Configure succeeded: ${INNER_BUILD_DIR}"
