#!/bin/bash
# \author Hans J. Johnson
#
# This script provides some notes for facilitating builds
# of BRAINSTools for specialized purposes.
#
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

BRAINS_DOC_DIR=$(dirname "${SCRIPT_DIR}")
BRAINS_SRC_DIR=$(dirname $(dirname "${SCRIPT_DIR}"))

cp CMakeUserPresets.json "${BRAINS_SRC_DIR}"

pushd "${BRAINS_SRC_DIR}" || exit

cmake -S "${BRAINS_SRC_DIR}" --preset brainstools_support
ninja -C "${BRAINS_SRC_DIR}/cmake-release-support"
cmake -S "${BRAINS_SRC_DIR}" --preset brainstools_release
ninja -C "${BRAINS_SRC_DIR}/cmake-release-build"

popd || exit
