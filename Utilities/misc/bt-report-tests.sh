#!/usr/bin/env bash
# bt-report-tests.sh — parse a CTest log and print a structured test report.
#
# Works with both serial and parallel (-j N) CTest output because it matches
# the result lines that CTest always emits:
#
#   N/Total Test  #ID: name ....   Passed  T sec
#   N/Total Test  #ID: name ....   Failed  T sec
#   N/Total Test  #ID: name ....***Not Run T sec
#
# Usage:
#   bt-report-tests.sh [OPTIONS] [LOGFILE]
#
# Arguments:
#   LOGFILE   CTest log to parse (default: TestResults/latest.log in CWD,
#             or stdin if "-" is passed)
#
# Options:
#   -n <N>    Number of slowest tests to list  (default: 20)
#   -o <FILE> Write tab-separated raw data to FILE (num\ttime\tstatus\tname)
#   -h        Show this help
#
# Exit code:
#   0   all tests passed (or no failures found in log)
#   1   one or more failures / not-run tests detected

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
TOP_N=20
RAW_OUT=""
LOGFILE=""

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -n) TOP_N="$2";   shift 2 ;;
    -o) RAW_OUT="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,/^set -/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    -*) echo "Unknown option: $1" >&2; exit 1 ;;
    *)  LOGFILE="$1"; shift ;;
  esac
done

if [[ -z "${LOGFILE}" ]]; then
  LOGFILE="${PWD}/TestResults/latest.log"
fi

if [[ "${LOGFILE}" == "-" ]]; then
  TMP_LOG=$(mktemp)
  cat > "${TMP_LOG}"
  LOGFILE="${TMP_LOG}"
  trap 'rm -f "${TMP_LOG}"' EXIT
fi

if [[ ! -f "${LOGFILE}" ]]; then
  echo "ERROR: log file not found: ${LOGFILE}" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Parse — emit tab-separated: num<TAB>time<TAB>status<TAB>name
# Handles both "Test #N:" (serial) and "Test  #N:" (parallel -j) spacing.
# ---------------------------------------------------------------------------
RAW=$(perl - "${LOGFILE}" <<'EOF'
  use strict; use warnings;
  open(my $fh, '<', $ARGV[0]) or die $!;
  while (<$fh>) {
    if (/Test\s+#(\d+):\s+(.+?)\s*\.+\s*(Passed|Failed|Not Run)\s+([\d.]+)\s+sec/) {
      printf "%d\t%.2f\t%s\t%s\n", $1, $4, $3, $2;
    }
  }
EOF
)

if [[ -z "${RAW}" ]]; then
  echo "ERROR: no test result lines found in: ${LOGFILE}" >&2
  exit 1
fi

# Optionally write raw data file
if [[ -n "${RAW_OUT}" ]]; then
  echo "${RAW}" | sort -k1 -n > "${RAW_OUT}"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
TOTAL=$(echo "${RAW}"    | wc -l | tr -d ' ')
PASSED=$(echo "${RAW}"   | awk -F'\t' '$3=="Passed"'   | wc -l | tr -d ' ')
FAILED=$(echo "${RAW}"   | awk -F'\t' '$3=="Failed"'   | wc -l | tr -d ' ')
NOTRUN=$(echo "${RAW}"   | awk -F'\t' '$3=="Not Run"'  | wc -l | tr -d ' ')
TOTAL_S=$(echo "${RAW}"  | awk -F'\t' '{s+=$2} END{printf "%.1f", s}')
WALL_S=$(grep -E "^Total Test time" "${LOGFILE}" 2>/dev/null | tail -1 \
         | grep -Eo '[0-9]+\.[0-9]+' || echo "n/a")

echo "════════════════════════════════════════════════════════════════"
echo " BRAINSTools Test Report"
echo " Log: ${LOGFILE}"
echo "════════════════════════════════════════════════════════════════"
printf " %-12s %d\n"  "Total:"    "${TOTAL}"
printf " %-12s %d\n"  "Passed:"   "${PASSED}"
printf " %-12s %d\n"  "Failed:"   "${FAILED}"
printf " %-12s %d\n"  "Not Run:"  "${NOTRUN}"
printf " %-12s %s s\n" "Sum(test):"  "${TOTAL_S}"
printf " %-12s %s s\n" "Wall clock:" "${WALL_S}"
echo "════════════════════════════════════════════════════════════════"

# ---------------------------------------------------------------------------
# Top-N slowest
# ---------------------------------------------------------------------------
echo ""
echo " Top ${TOP_N} slowest tests:"
echo " ┌──────┬──────────────┬────────────┬──────────────────────────────────────────────┐"
printf " │ %4s │ %12s │ %-10s │ %-44s │\n" "  #" "    Time (s)" "Status" "Test name"
echo " ├──────┼──────────────┼────────────┼──────────────────────────────────────────────┤"
echo "${RAW}" | sort -t$'\t' -k2 -rn | head -"${TOP_N}" | \
  awk -F'\t' '{printf " │ %4d │ %12.2f │ %-10s │ %-44.44s │\n", $1, $2, $3, $4}'
echo " └──────┴──────────────┴────────────┴──────────────────────────────────────────────┘"

# ---------------------------------------------------------------------------
# Failures / Not Run
# ---------------------------------------------------------------------------
FAIL_LINES=$(echo "${RAW}" | awk -F'\t' '$3 != "Passed"')
echo ""
if [[ -z "${FAIL_LINES}" ]]; then
  echo " ✅  All ${TOTAL} tests passed."
  EXIT_CODE=0
else
  echo " ❌  Failures / Not Run:"
  echo " ┌──────┬──────────────┬────────────┬──────────────────────────────────────────────┐"
  printf " │ %4s │ %12s │ %-10s │ %-44s │\n" "  #" "    Time (s)" "Status" "Test name"
  echo " ├──────┼──────────────┼────────────┼──────────────────────────────────────────────┤"
  echo "${FAIL_LINES}" | sort -k1 -n | \
    awk -F'\t' '{printf " │ %4d │ %12.2f │ %-10s │ %-44.44s │\n", $1, $2, $3, $4}'
  echo " └──────┴──────────────┴────────────┴──────────────────────────────────────────────┘"
  EXIT_CODE=1
fi

echo ""
exit "${EXIT_CODE}"
