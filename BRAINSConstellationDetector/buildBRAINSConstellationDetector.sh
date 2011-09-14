#!/bin/bash

BUILDNAME=$(uname)_$(uname -m)-$(hostname -s)

SOURCE_DIR=$(dirname $0)
if [ "${SOURCE_DIR}" == "." ]; then
  SOURCE_DIR=$(pwd)
fi

if [ $# -lt 1 ]; then
  ABI="FAST"
else
  ABI=$1
fi
ITK_BUILD_NAME=$(uname).$(uname -m).${ABI}

## Valid types are Experimental, Continous, Nightly
if [ $# -lt 2 ]; then
  BUILDTYPE=Experimental
else
  BUILDTYPE=$2
fi

CC=/usr/bin/gcc-4.2
CXX=/usr/bin/g++-4.2
if [ ! -f ${CC} ] || [ ! -f ${CXX} ]; then
CC=gcc
CXX=g++
fi
case ${ABI} in
  "PROFILE")
    CFLAGS="-Wall -Wstrict-prototypes -fprofile-arcs -ftest-coverage -pg  -UNDEBUG"
    CXXFLAGS="-Wall  -fprofile-arcs -ftest-coverage -pg -UNDEBUG"
    ;;
  "DEBUG")
    CFLAGS="-Wall -Wstrict-prototypes -g"
    CXXFLAGS="-Wall  -g"
    ;;
  "FAST")
    CFLAGS="-DNDEBUG -O3 -msse -mmmx -msse2 -msse3"
    CXXFLAGS="-DNDEBUG -O3 -msse -mmmx -msse2 -msse3"
    ;;
  *)
    echo "INVALID ABI GIVEN"
    exit -1;
esac

COMPILE_DIR=$(dirname ${SOURCE_DIR})/ART-COMPILE/
ABI_DIR=${COMPILE_DIR}/$(uname)_$(uname -m)-${ABI}
mkdir -p ${ABI_DIR}

if [ 0 -eq 1 ]; then
  #################################################################################
  #Get and build Insight-3.8.0
  ITK_SOURCE=${COMPILE_DIR}/InsightToolkit-3.8.0
  ITK_BUILD=${ABI_DIR}/InsightToolkit-3.8.0-build
  pushd ${COMPILE_DIR}
  if [ ! -d ${ITK_SOURCE} ]; then
    WGETBIN=$(which wget)
    if [ "x${WGETBIN}" == "x" ] || [ ! -f ${WGETBIN} ]; then
      echo "ERROR: Can not build without wget to retrive ITK."
      exit -1;
    fi
    wget http://voxel.dl.sourceforge.net/sourceforge/itk/InsightToolkit-3.8.0.tar.gz
    tar -xzvf InsightToolkit-3.8.0.tar.gz
    popd
  fi
  mkdir -p ${ITK_BUILD}
  pushd ${ITK_BUILD}
  ##NOTE:  Using cmake and all comand line options.  Normally ccmake would be used.
  CC=${CC} CXX=${CXX} CFLAGS=${CFLAGS} CXXFLAGS=${CXXFLAGS} cmake ${ITK_SOURCE} -DBUILD_EXAMPLES:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=OFF -DBUILDNAME:STRING=${ITK_BUILD_NAME} -DCOVERAGE_COMMAND:FILEPATH=/usr/bin/gcov-4.2
  if [ $? -ne 0 ]; then
    echo "ERROR in building ITK"
    exit -1
  fi
else
  #################################################################################
  #Get and build InsightToolkit-CVS
  LOCAL_PATH=$(dirname $0)
  echo "${LOCAL_PATH}"
  if [ "${LOCAL_PATH}" == "." ]; then
    LOCAL_PATH=$(pwd)
  fi
  echo "${LOCAL_PATH}"
  ITK_TARBALL=${LOCAL_PATH}/InsightToolkit-CVS.tar.gz
  ITK_SOURCE=${COMPILE_DIR}/Insight
  ITK_BUILD=${ABI_DIR}/InsightToolkit-CVS-build
  if [ ! -f ${ITK_SOURCE}/CMakeLists.txt ] || [ ${LOCAL_PATH}/InsightToolkit-CVS.tar.gz -nt ${ITK_SOURCE}/CMakeLists.txt ]; then
    mkdir -p ${ITK_SOURCE}
    pushd ${COMPILE_DIR}
    tar -xvzf ${ITK_TARBALL}
    touch ${ITK_SOURCE}/CMakeLists.txt
    popd
  fi
  mkdir -p ${ITK_BUILD}
  pushd ${ITK_BUILD}
  ##NOTE:  Using cmake and all comand line options.  Normally ccmake would be used.
  CC=${CC} CXX=${CXX} CFLAGS=${CFLAGS} CXXFLAGS=${CXXFLAGS} cmake ${ITK_SOURCE} \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=OFF \
    -DBUILDNAME:STRING=${ITK_BUILD_NAME} \
    -DCOVERAGE_COMMAND:FILEPATH=/usr/bin/gcov-4.2 \
    -DITK_USE_ORIENTED_IMAGE_DIRECTION:BOOL=ON \
    -DITK_USE_REVIEW:BOOL=ON \
    -DITK_USE_TRANSFORM_IO_FACTORIES:BOOL=ON \
    -DITK_USE_OPTIMIZED_REGISTRATION_METHODS:BOOL=ON \
    -DITK_USE_ORIENTED_IMAGE_DIRECTION:BOOL=ON \
    -DITK_IMAGE_BEHAVES_AS_ORIENTED_IMAGE:BOOL=ON
  if [ $? -ne 0 ]; then
    echo "ERROR in configuring ITK"
    exit -1
  fi
fi
## Auto-determine the number of processors on this computer
case "$(uname)" in
  "Linux")
  maxproc=$(grep processor /proc/cpuinfo | wc | awk '{print $1}')
  ;;
  "Darwin")
  maxproc=$(sysctl -n hw.ncpu)
  ;;
  *)
  echo "Platform not recognized"
  maxproc=1;
esac
if [ "x${maxproc}" == "x" ] || [ ${maxproc} -lt 1 ]; then
  maxproc=1;
fi
make -j${maxproc}
popd

#################################################################################
#Build ART
mkdir -p ${ABI_DIR}/ART
pushd ${ABI_DIR}/ART
##NOTE:  Using cmake and all comand line options.  Normally ccmake would be used.
CC=${CC} CXX=${CXX} CFLAGS=${CFLAGS} CXXFLAGS=${CXXFLAGS} cmake ${SOURCE_DIR}/ -DBUILD_EXAMPLES:BOOL=ON -DBUILD_TESTING:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=OFF -DITK_DIR:PATH=${ITK_BUILD} -DABI:STRING=${ABI} -DBUILDNAME:STRING=${ITK_BUILD_NAME}  -DCOVERAGE_COMMAND:FILEPATH=/usr/bin/gcov-4.2


case "$(uname)" in
  "Linux")
  maxproc=$(grep processor /proc/cpuinfo | wc | awk '{print $1}')
  ;;
  "Darwin")
  maxproc=$(sysctl -n hw.ncpu)
  ;;
  *)
  echo "Platform not recognized"
  maxproc=1;
esac
make ${BUILDTYPE}
popd


