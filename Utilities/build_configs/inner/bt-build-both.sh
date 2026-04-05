#!/usr/bin/env bash
# bt-build-both.sh
# Configure (if needed) and build both Debug and Release BRAINSTools inner
# builds, reusing the ExternalProject trees from an existing SuperBuild.
#
# Usage:
#   bt-build-both.sh [OPTIONS]
#
# Options:
#   -b, --build-dir DIR    SuperBuild binary directory   (default: auto-detect)
#   -s, --source-dir DIR   BRAINSTools source root       (default: auto-detect)
#   -j, --jobs N           Parallel jobs for ninja       (default: nproc/2)
#   -f, --force            Re-configure even if already configured
#   -c, --configure-only   Configure but do not build
#   -d, --debug-only       Build only the Debug inner build
#   -r, --release-only     Build only the Release inner build
#   -h, --help             Show this message
#
# Environment:
#   BRAINSTOOLS_BUILD_DIR  overrides --build-dir auto-detection
#   BRAINSTOOLS_SOURCE_DIR overrides --source-dir auto-detection

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIGURE_SCRIPT="${SCRIPT_DIR}/bt-configure-inner.sh"

# -----------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------
FORCE_FLAG=""
CONFIGURE_ONLY=0
BUILD_DEBUG=1
BUILD_RELEASE=1
EXTRA_ARGS=()

# Parallel jobs: default to half the logical CPUs (inner builds are CPU-heavy)
if command -v nproc &>/dev/null; then
  DEFAULT_JOBS=$(( $(nproc) / 2 ))
else
  DEFAULT_JOBS=$(( $(sysctl -n hw.logicalcpu 2>/dev/null || echo 4) / 2 ))
fi
JOBS="${DEFAULT_JOBS}"

# -----------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--build-dir)    EXTRA_ARGS+=(--build-dir "$2"); shift 2 ;;
    -s|--source-dir)   EXTRA_ARGS+=(--source-dir "$2"); shift 2 ;;
    -j|--jobs)         JOBS="$2"; shift 2 ;;
    -f|--force)        FORCE_FLAG="--force"; shift ;;
    -c|--configure-only) CONFIGURE_ONLY=1; shift ;;
    -d|--debug-only)   BUILD_RELEASE=0; shift ;;
    -r|--release-only) BUILD_DEBUG=0; shift ;;
    -h|--help)
      sed -n '2,/^set -/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# -----------------------------------------------------------------------
# Helper: resolve inner build dir from superbuild cache
# -----------------------------------------------------------------------
_inner_build_dir() {
  local build_type="$1"
  # Re-use the same auto-detect logic: find SUPERBUILD_DIR from EXTRA_ARGS or env
  local sb_dir="${BRAINSTOOLS_BUILD_DIR:-}"
  for i in "${!EXTRA_ARGS[@]}"; do
    if [[ "${EXTRA_ARGS[$i]}" == "--build-dir" ]]; then
      sb_dir="${EXTRA_ARGS[$((i+1))]}"
    fi
  done
  if [[ -z "${sb_dir}" ]]; then
    local auto_source="${BRAINSTOOLS_SOURCE_DIR:-$(cd "${SCRIPT_DIR}/../../.." && pwd)}"
    sb_dir="${auto_source}/build"
  fi
  echo "${sb_dir}/BRAINSTools-${build_type}-EPRelease-build"
}

# -----------------------------------------------------------------------
# Configure
# -----------------------------------------------------------------------
if [[ "${BUILD_DEBUG}" -eq 1 ]]; then
  echo ">>> Configure: Debug"
  "${CONFIGURE_SCRIPT}" --type Debug "${EXTRA_ARGS[@]}" ${FORCE_FLAG}
fi

if [[ "${BUILD_RELEASE}" -eq 1 ]]; then
  echo ">>> Configure: Release"
  "${CONFIGURE_SCRIPT}" --type Release "${EXTRA_ARGS[@]}" ${FORCE_FLAG}
fi

[[ "${CONFIGURE_ONLY}" -eq 1 ]] && { echo "Configure-only mode — done."; exit 0; }

# -----------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------
build_one() {
  local build_type="$1"
  local dir
  dir="$(_inner_build_dir "${build_type}")"

  echo ""
  echo ">>> Build: ${build_type}  (${dir})"
  echo "    Jobs: ${JOBS}"

  local start_s
  start_s=$(date +%s)

  ninja -C "${dir}" -j "${JOBS}"

  local elapsed=$(( $(date +%s) - start_s ))
  echo "    Build ${build_type} complete in ${elapsed}s"
}

BUILD_FAILED=0

if [[ "${BUILD_DEBUG}" -eq 1 ]]; then
  build_one Debug || { echo "ERROR: Debug build failed" >&2; BUILD_FAILED=1; }
fi

if [[ "${BUILD_RELEASE}" -eq 1 ]]; then
  build_one Release || { echo "ERROR: Release build failed" >&2; BUILD_FAILED=1; }
fi

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
echo ""
echo "=== Build summary ==="
if [[ "${BUILD_DEBUG}" -eq 1 ]]; then
  dir="$(_inner_build_dir Debug)"
  bin_count=$(find "${dir}" -maxdepth 4 -type f -perm -u+x ! -path "*/CMakeFiles/*" ! -name "*.cmake" ! -name "*.sh" 2>/dev/null | wc -l | tr -d ' ')
  echo "  Debug   : ${dir}"
  echo "            ${bin_count} executable targets"
fi
if [[ "${BUILD_RELEASE}" -eq 1 ]]; then
  dir="$(_inner_build_dir Release)"
  bin_count=$(find "${dir}" -maxdepth 4 -type f -perm -u+x ! -path "*/CMakeFiles/*" ! -name "*.cmake" ! -name "*.sh" 2>/dev/null | wc -l | tr -d ' ')
  echo "  Release : ${dir}"
  echo "            ${bin_count} executable targets"
fi

exit "${BUILD_FAILED}"
