#!/bin/bash
# \author Hans J. Johnson
#
# This is a wrapper around the non-standard, often confusing
# tbb build environment
#
# It unfortunately requires make, python, and cmake in a
# wierd configuration to get a standard build and install
# layout.
#
# Usage:  buld_tbb.sh -b build_directory -p install_prefix
#    EG:  ./build_tbb.sh -b /tmp -p /home/myaccout/local
set -ev

while getopts ":b:p:" opt; do
  case ${opt} in
    b ) # build
      BUILD_CACHE=$OPTARG
      ;;
    p ) # install
      INSTALL_PREFIX=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-h] [-t]"
      ;;
  esac
done

if [ -z ${BUILD_CACHE+x} ]; then
    echo "Usage was wrong, -b is a required option"
    exit 255
fi
if [ -z ${INSTALL_PREFIX+x} ]; then
    INSTALL_PREFIX=${BUILD_CACHE}/tbb-install
fi

SRC_DIR=${BUILD_CACHE}

#SYSTEM_NAME=Darwin
TBB_VERSION_FILE=${SRC_DIR}/include/tbb/tbb_stddef.h

mkdir -p "${BUILD_CACHE}"
mkdir -p "${INSTALL_PREFIX}"

pushd "${BUILD_CACHE}"

python3 "${SRC_DIR}/build/build.py"  \
       --tbbroot "${SRC_DIR}" \
       --prefix "${INSTALL_PREFIX}" \
       --install-libs --install-devel --install-docs \
  > /tmp/tbb_logger 2>&1

#cat /tmp/tbb_logger

cmake -DINSTALL_DIR="${INSTALL_PREFIX}/lib/cmake/tbb" \
      -DLIB_REL_PATH=../../../lib \
      -DINC_REL_PATH=../../../include \
      -DSYSTEM_NAME="$(uname)" \
      -DTBB_VERSION_FILE="${TBB_VERSION_FILE}" \
      -P "${SRC_DIR}/cmake/tbb_config_installer.cmake"


DO_TESTING=1
if [[ ${DO_TESTING} -eq 1 ]]; then

TEST_DIR=${BUILD_CACHE}/tbb_test
mkdir -p "${TEST_DIR}"

cat > "${TEST_DIR}/CMakeLists.txt" << EOF
# Put this cmake file in a directory, run "cmake ."  finding TBB results in a FATAL_ERROR from cmake
cmake_minimum_required(VERSION 3.0)
# tbb;tbbmalloc;tbbmalloc_proxy
find_package(TBB REQUIRED COMPONENTS tbbmalloc CONFIG)       #Find one component - test
find_package(TBB REQUIRED COMPONENTS tbbmalloc_proxy CONFIG) #Find one component - test
find_package(TBB REQUIRED COMPONENTS tbb CONFIG)             #Find one component - test
find_package(TBB REQUIRED CONFIG) # Should be able to find_package many times
if(NOT TBB_FOUND)
  message(FATAL_ERROR "NO TBB")
else()
  #message(STATUS "COMPONENTS FOUND: TBB_IMPORTED_TARGETS=:\${TBB_IMPORTED_TARGETS}: TBB_tbb_FOUND=:\${TBB_tbb_FOUND}: TBB_tbbmalloc_FOUND=:\${TBB_tbbmalloc_FOUND}: TBB_tbbmalloc_proxy_FOUND=:\${TBB_tbbmalloc_proxy_FOUND}:")
  message(STATUS "FOUND TBB")
endif()
add_executable(tbb_test test.cpp)
target_link_libraries(tbb_test TBB::tbb)
EOF

cat > "${TEST_DIR}/test.cpp" << EOF
#include <tbb/task_scheduler_init.h>
#include <iostream>
int main()
{
  std::cout << "Default number of threads: " << tbb::task_scheduler_init::default_num_threads() << std::endl;;
  return 0;
}
EOF

mkdir -p "${TEST_DIR}-bld"
pushd "${TEST_DIR}-bld"
cmake -DCMAKE_BUILD_TYPE:STRING=Release -DTBB_DIR:PATH="${INSTALL_PREFIX}/lib/cmake/tbb"  "${TEST_DIR}"
make

./tbb_test
tbb_test_status=$?
if [[ "${tbb_test_status}"  -eq 0 ]]; then
  echo "SUCCESSFUL TESTING"
  exit 0
else
  echo "FAILED TESTING"
  exit 255
fi

fi
