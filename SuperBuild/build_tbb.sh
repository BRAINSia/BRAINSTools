#!/bin/bash -eax
# \author Hans J. Johnson
#
# This is a wrapper around the non-standard, often confusing
# tbb build environment
#
# It unfortunately requires make, python, and cmake in a
# wierd configuration to get a standard build and install
# layout.
#
# Usage:  build_tbb.sh -b build_directory -p install_prefix
#    EG:  ./build_tbb.sh -b /tmp -p /home/myaccout/local
set -ev

while getopts ":c:x:i:b:p:d:" opt; do
  case ${opt} in
    c ) # compiler_id
      CC=$OPTARG
      ;;
    x ) # compiler_id
      CXX=$OPTARG
      ;;
    i ) # compiler_id
      TBB_COMPILERID=$OPTARG
      ;;
    b ) # build
      BUILD_CACHE=$OPTARG
      ;;
    p ) # install
      INSTALL_PREFIX=$OPTARG
      ;;
    d ) # cmake -D flags
      CMAKE_CMD_FLAGS="$OPTARG"
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

export CXX=${CXX}
export CC=${CC}
if [ "x$(uname -s)" = "xDarwin" ]; then
  # Building TBB on mac is incredibly frustrating!!
  export CFLAGS="-isysroot $(xcrun --sdk macosx --show-sdk-path) ${CFLAGS}"
  export CXXFLAGS="-isysroot $(xcrun --sdk macosx --show-sdk-path) ${CXXFLAGS}"
fi
CFLAGS=${CFLAGS} CXXFLAGS=${CXXFLAGS} CXX=${CXX} CC=${CC} python3 "${SRC_DIR}/build/build.py"  \
       --tbbroot "${SRC_DIR}" \
       --prefix "${INSTALL_PREFIX}" \
       --install-libs --install-devel --install-docs \
       --build-args compiler=${TBB_COMPILERID} \
  > /tmp/tbb_logger 2>&1
python_build_status=$?
if [ $python_build_status -ne 0 ]; then
  cat /tmp/tbb_logger
  exit ${python_build_status}
fi

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

project(TestFindingTBB
       VERSION 0.0.0.1
       LANGUAGES CXX C
)

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
//#include <tbb/task_scheduler_init.h> //deprecated
#include "tbb/task_arena.h"
#include <iostream>
int main()
{
  //std::cout << "Default number of threads: " << tbb::task_scheduler_init::default_num_threads() << std::endl;;
  tbb::task_arena default_arena;
  std::cout << " " <<  default_arena.max_concurrency() << std::endl;
  return 0;
}
EOF

  mkdir -p "${TEST_DIR}-bld"
  pushd "${TEST_DIR}-bld"
  echo "====="
  cmake -DCMAKE_BUILD_TYPE:STRING=Release \
    -DTBB_DIR:PATH=${INSTALL_PREFIX}/lib/cmake/tbb  \
     ${CMAKE_CMD_FLAGS} \
    "${TEST_DIR}"
  if [[ -f Makefile ]];then
    echo "Attempting to build with a Makefile"
    make
  else
    echo "Attempting to build with a ninja"
    ninja
  fi
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
