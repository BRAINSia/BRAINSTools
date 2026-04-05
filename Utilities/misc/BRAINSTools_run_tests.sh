#!/usr/bin/env bash
# Run BRAINSTools CTests with controlled parallelism.
#
# Tests declare PROCESSORS 3 in CMakeLists so CTest manages scheduling
# natively.  The default parallel job count is:
#
#   parallel_jobs = hardware_cpus / 3
#
# Results are saved to TestResults/ in the current build directory as a
# timestamped JUnit XML file (requires CMake >= 3.21) and a plain-text log.
# A symlink TestResults/latest.junit.xml always points to the most recent run.
# To address failing tests without re-running everything:
#
#   ../../Utilities/misc/BRAINSTools_run_tests.sh -- --rerun-failed
#
# Usage:
#   cd <build-dir>
#   ../../Utilities/misc/BRAINSTools_run_tests.sh [options]
#
# Options:
#   -j <N>   parallel ctest jobs   (default: hw_cpus / 3)
#   -r <R>   ctest -R regex filter (default: none)
#   --       pass remaining args directly to ctest
#
# Examples:
#   ../../Utilities/misc/BRAINSTools_run_tests.sh
#   ../../Utilities/misc/BRAINSTools_run_tests.sh -j 4
#   ../../Utilities/misc/BRAINSTools_run_tests.sh -r BRAINSFit
#   ../../Utilities/misc/BRAINSTools_run_tests.sh -- --rerun-failed

set -euo pipefail

# ---------------------------------------------------------------------------
# Detect hardware CPU count (Linux and macOS)
# ---------------------------------------------------------------------------
if command -v nproc &>/dev/null; then
    HW_CPUS=$(nproc)
elif command -v sysctl &>/dev/null; then
    HW_CPUS=$(sysctl -n hw.ncpu)
else
    HW_CPUS=4
    echo "WARNING: could not detect CPU count, defaulting to ${HW_CPUS}" >&2
fi

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
PROCESSORS_PER_TEST=3
PARALLEL_JOBS=$(( HW_CPUS / PROCESSORS_PER_TEST ))
REGEX=""
EXTRA_ARGS=()

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j)  PARALLEL_JOBS="$2"; shift 2 ;;
        -r)  REGEX="$2"; shift 2 ;;
        --)  shift; EXTRA_ARGS+=("$@"); break ;;
        *)   echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ "${PARALLEL_JOBS}" -lt 1 ]]; then
    PARALLEL_JOBS=1
fi

# ---------------------------------------------------------------------------
# Set up timestamped results directory
# ---------------------------------------------------------------------------
TIMESTAMP=$(date +%Y%m%dT%H%M%S)
RESULTS_DIR="${PWD}/TestResults"
mkdir -p "${RESULTS_DIR}"
JUNIT_FILE="${RESULTS_DIR}/${TIMESTAMP}.junit.xml"
LOG_FILE="${RESULTS_DIR}/${TIMESTAMP}.log"

# ---------------------------------------------------------------------------
# Build ctest command
# ---------------------------------------------------------------------------
CTEST_CMD=(
    ctest
    -j "${PARALLEL_JOBS}"
    --output-on-failure
    --output-junit "${JUNIT_FILE}"
)
[[ -n "${REGEX}" ]] && CTEST_CMD+=(-R "${REGEX}")
CTEST_CMD+=("${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}")

# ---------------------------------------------------------------------------
# Report and run
# ---------------------------------------------------------------------------
echo "Hardware CPUs         : ${HW_CPUS}"
echo "PROCESSORS per test   : ${PROCESSORS_PER_TEST}"
echo "Parallel ctest jobs   : ${PARALLEL_JOBS}"
[[ -n "${REGEX}" ]] && echo "Test regex filter     : ${REGEX}"
echo "Working directory     : $(pwd)"
echo "JUnit XML             : ${JUNIT_FILE}"
echo "Text log              : ${LOG_FILE}"
echo "Command               : ${CTEST_CMD[*]}"
echo ""

# Run ctest, capturing output to log while also displaying it
CTEST_EXIT=0
"${CTEST_CMD[@]}" 2>&1 | tee "${LOG_FILE}" || CTEST_EXIT=${PIPESTATUS[0]}

# Update the stable symlink to the latest run
ln -sf "${TIMESTAMP}.junit.xml" "${RESULTS_DIR}/latest.junit.xml"
ln -sf "${TIMESTAMP}.log"       "${RESULTS_DIR}/latest.log"

echo ""
echo "Results saved:"
echo "  XML : ${JUNIT_FILE}"
echo "  log : ${LOG_FILE}"
echo "  symlinks: ${RESULTS_DIR}/latest.junit.xml  ${RESULTS_DIR}/latest.log"

# ---------------------------------------------------------------------------
# Print structured report (summary + slowest + failures)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo ""
"${SCRIPT_DIR}/bt-report-tests.sh" "${LOG_FILE}" || true

exit "${CTEST_EXIT}"
