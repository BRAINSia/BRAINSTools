#!/usr/bin/env bash
# bt-profile-tests.sh
# Profile BRAINSTools tests using macOS `sample`.
#
# Usage:
#   bt-profile-tests.sh [OPTIONS]
#
# Options:
#   -b, --build-dir DIR   Inner build directory
#                         (default: auto-detect Profiling build, fall back to Release)
#   -o, --output-dir DIR  Directory for sample output files  (default: /tmp/bt-profiles)
#   -s, --stub NAME       Run only the named stub (e.g. BRAINSABCSmallTest).
#                         May be specified multiple times.  Default: all stubs.
#   -l, --list            List available stubs and exit.
#   -h, --help            Show this message
#
# Stub scripts live in:
#   Utilities/misc/bt-profile-stubs/
#
# Each stub is a self-contained bash file that defines:
#   _stub_name  – test name used for the sample file
#   _stub_cmd   – bash array of the full command to profile
#   _stub_out   – path to the .sample.txt output file
#
# To profile a single test interactively, run the stub directly:
#   bt-profile-stubs/brainsabc-small.sh --build-dir DIR --data-dir DIR
#
# Requirements:
#   macOS only — uses the `sample` CLI tool bundled with Xcode Command Line Tools.
#   Build must have been made with BUILD_PROFILING=ON for accurate symbol resolution.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUBS_DIR="${SCRIPT_DIR}/bt-profile-stubs"
SOURCE_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SUPERBUILD_DIR="${BRAINSTOOLS_BUILD_DIR:-${SOURCE_DIR}/build}"
OUTPUT_DIR="/tmp/bt-profiles"

# -----------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------
INNER_BUILD_DIR=""
SELECTED_STUBS=()
LIST_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--build-dir)  INNER_BUILD_DIR="$2"; shift 2 ;;
    -o|--output-dir) OUTPUT_DIR="$2";      shift 2 ;;
    -s|--stub)       SELECTED_STUBS+=("$2"); shift 2 ;;
    -l|--list)       LIST_ONLY=1; shift ;;
    -h|--help)
      sed -n '2,/^set -/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# -----------------------------------------------------------------------
# Discover stubs
# -----------------------------------------------------------------------
mapfile -t ALL_STUB_FILES < <(find "${STUBS_DIR}" -name "*.sh" | sort)

if [[ ${LIST_ONLY} -eq 1 ]]; then
  echo "Available stubs (in ${STUBS_DIR}):"
  for f in "${ALL_STUB_FILES[@]}"; do
    # Source stub in a subshell to read _stub_name without side-effects
    _stub_name=""
    # shellcheck source=/dev/null
    ( source "$f" 2>/dev/null; echo "  ${_stub_name}  ← $(basename "$f")" ) || true
  done
  exit 0
fi

# -----------------------------------------------------------------------
# Auto-detect inner build directory
# -----------------------------------------------------------------------
if [[ -z "${INNER_BUILD_DIR}" ]]; then
  for candidate in \
    "${SUPERBUILD_DIR}/BRAINSTools-Release-Profiling-EPRelease-build" \
    "${SUPERBUILD_DIR}/BRAINSTools-Release-EPRelease-build"
  do
    if [[ -f "${candidate}/BRAINSFit/TestSuite/BRAINSFitTestDriver" ]]; then
      INNER_BUILD_DIR="${candidate}"
      break
    fi
  done
fi

if [[ -z "${INNER_BUILD_DIR}" || ! -f "${INNER_BUILD_DIR}/BRAINSFit/TestSuite/BRAINSFitTestDriver" ]]; then
  echo "ERROR: Could not find an inner build directory with test drivers. Use --build-dir or build first." >&2
  exit 1
fi

# Shared ExternalData directory (always from the Release build)
DATA_DIR="${SUPERBUILD_DIR}/BRAINSTools-Release-EPRelease-build/ExternalData/TestData"
if [[ ! -d "${DATA_DIR}" ]]; then
  echo "ERROR: ExternalData not found at ${DATA_DIR}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "=== BRAINSTools Performance Profiler ==="
echo "  Inner build : ${INNER_BUILD_DIR}"
echo "  Data dir    : ${DATA_DIR}"
echo "  Output dir  : ${OUTPUT_DIR}"
echo ""

# -----------------------------------------------------------------------
# Helper: run one command under `sample`, collect output
# -----------------------------------------------------------------------
_profile_test() {
  local test_name="$1"
  local sample_file="${OUTPUT_DIR}/${test_name}.sample.txt"
  shift
  local -a cmd=("$@")

  echo "--- Profiling: ${test_name} ---"
  echo "  Command: ${cmd[*]}"
  echo ""

  "${cmd[@]}" &
  local pid=$!

  sample "${pid}" 600 -wait -f "${sample_file}" 2>/dev/null || true

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
# Source stubs and run
# -----------------------------------------------------------------------
RAN=0
for stub_file in "${ALL_STUB_FILES[@]}"; do
  # Load stub: defines _stub_name, _stub_cmd[], _stub_out
  _stub_name="" _stub_out="" _stub_cmd=()
  # shellcheck source=/dev/null
  source "${stub_file}"

  # Filter by --stub selection if given
  if [[ ${#SELECTED_STUBS[@]} -gt 0 ]]; then
    match=0
    for sel in "${SELECTED_STUBS[@]}"; do
      [[ "${_stub_name}" == "${sel}" ]] && match=1 && break
    done
    [[ ${match} -eq 0 ]] && continue
  fi

  _profile_test "${_stub_name}" "${_stub_cmd[@]}"
  (( RAN++ )) || true
done

if [[ ${RAN} -eq 0 ]]; then
  echo "ERROR: No matching stubs found." >&2
  exit 1
fi

# -----------------------------------------------------------------------
# Hotspot report for every sample file in OUTPUT_DIR
# -----------------------------------------------------------------------
echo ""
echo "================================================================"
echo " HOTSPOT REPORT"
echo "================================================================"
echo ""

for stub_file in "${ALL_STUB_FILES[@]}"; do
  _stub_name="" _stub_out=""
  # shellcheck source=/dev/null
  source "${stub_file}" 2>/dev/null || true
  [[ -z "${_stub_name}" ]] && continue

  if [[ ${#SELECTED_STUBS[@]} -gt 0 ]]; then
    match=0
    for sel in "${SELECTED_STUBS[@]}"; do
      [[ "${_stub_name}" == "${sel}" ]] && match=1 && break
    done
    [[ ${match} -eq 0 ]] && continue
  fi

  _top_hotspots "${_stub_name}" "${OUTPUT_DIR}/${_stub_name}.sample.txt"
done

echo "Sample files saved to: ${OUTPUT_DIR}/"
echo ""
echo "To explore interactively, open in Instruments:"
echo "  xcrun xctrace import --input <file> --output out.trace && open out.trace"
