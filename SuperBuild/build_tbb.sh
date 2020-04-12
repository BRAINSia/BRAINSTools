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
  echo "BUILDING FOR Darwin"
  export CFLAGS="-std=c11 -isysroot $(xcrun --sdk macosx --show-sdk-path) ${CFLAGS} -DTBB_DEPRECATED=0 -DTBB_USE_CAPTURED_EXCEPTION=0"
  export CXXFLAGS="-std=c++11 -isysroot $(xcrun --sdk macosx --show-sdk-path) ${CXXFLAGS} -DTBB_DEPRECATED=0 -DTBB_USE_CAPTURED_EXCEPTION=0"
else
  echo "BUILDING FOR LINUX"
  export CFLAGS="-std=c11 ${CFLAGS} -DTBB_DEPRECATED=0 -DTBB_USE_CAPTURED_EXCEPTION=0"
  export CXXFLAGS="-std=c++11 ${CXXFLAGS} -DTBB_DEPRECATED=0 -DTBB_USE_CAPTURED_EXCEPTION=0"
fi

python3 "${SRC_DIR}/build/build.py"  \
       --tbbroot "${SRC_DIR}" \
       --prefix "${INSTALL_PREFIX}" \
       --install-libs --install-devel --install-docs \
       --build-args "compiler='${TBB_COMPILERID}' CFLAGS='${CFLAGS}' CXXFLAGS='${CXXFLAGS}' CXX='${CXX}' CC='${CC}'" \
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
#include "tbb/task_arena.h"
#include <iostream>

#include "tbb/global_control.h"
// From tbb/examples/common/utility/get_default_num_threads.h
namespace tbb_utility
{
inline int
get_default_num_threads()
{
#if __TBB_MIC_OFFLOAD
#  pragma offload target(mic) out(default_num_threads)
#endif // __TBB_MIC_OFFLOAD
  static const size_t default_num_threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
  return static_cast<int>(default_num_threads);
}
} // namespace tbb_utility

int main()
{
  tbb::task_arena default_arena;
  std::cout << "Arena max concurrency: " <<  default_arena.max_concurrency() << std::endl;

  const int def_num_threads = tbb_utility::get_default_num_threads();
  std::cout << "Default number of threads: " << def_num_threads << std::endl;;

  const int initial_num_threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
    std::cout << "initial_num_threads: " << initial_num_threads << std::endl;
    if( initial_num_threads != def_num_threads )
    {
      std::cerr << "ERROR: initial_num_threads != def_num_threads " << std::endl;
      exit(100);
    }

  if( def_num_threads >=2 )
  {
    // Construct TBB static context with only 2 threads
    // https://software.intel.com/en-us/node/589744
    // Total parallelism that TBB can utilize
    // is limited by the current active global_control object
    // for the dynamic extension of the given scope.
    // ( instantiation of "global_control" object pushes the
    //   value onto active head of the FIFO stack for 'parameter'
    //   type, and destruction of the "global_control" object pops
    //   the 'parameter' type and returns the active value
    //   to it's original state.

    tbb::global_control test_tbb_global_context(
      tbb::global_control::max_allowed_parallelism,
      std::min<int>(tbb_utility::get_default_num_threads(), 2));
    const int local_num_threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
    std::cout << "Local num threads: " << local_num_threads << std::endl;
    if( local_num_threads != 2 )
    {
      std::cerr << "ERROR: local_num_threads != 2 " << std::endl;
      exit(101);
    }
  }
  const int restored_num_threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
  std::cout << "Restored num threads: " << restored_num_threads << std::endl;
  if( restored_num_threads != initial_num_threads )
    {
      std::cerr << "ERROR: restored_num_threads != initial_num_threads " << std::endl;
      exit(102);
    }
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
