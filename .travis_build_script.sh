#!/usr/bin/env bash

# -e: fail on error
# -v: show commands
# -x: show expanded commands
set -vex


SOURCE_DIR=$(pwd)
BUILD_DIR=${SOURCE_DIR}/cmake-travis-bld

mkdir ${BUILD_DIR}
cd ${BUILD_DIR}

cmake \
    -DCMAKE_CXX_STANDARD:STRING=11 \
    -DUSE_ANTS=BOOL=OFF \
    -DBRAINSTools_REQUIRES_VTK=OFF \
    -DUSE_BRAINSABC:BOOL=OFF \
    -DUSE_BRAINSDWICleanup:BOOL=OFF \
    -DUSE_BRAINSInitializedControlPoints:BOOL=OFF \
    -DUSE_BRAINSLabelStats:BOOL=OFF \
    -DUSE_BRAINSLandmarkInitializer:BOOL=OFF \
    -DUSE_BRAINSROIAuto:BOOL=OFF \
    -DUSE_BRAINSRefacer:BOOL=OFF \
    -DUSE_BRAINSResample:BOOL=OFF \
    -DUSE_BRAINSSnapShotWriter:BOOL=OFF \
    -DUSE_BRAINSStripRotation:BOOL=OFF \
    -DUSE_BRAINSTransformConvert:BOOL=OFF \
    -DUSE_ConvertBetweenFileFormats:BOOL=OFF \
    -DUSE_ImageCalculator:BOOL=OFF \
    -DUSE_ReferenceAtlas:BOOL=OFF \
    ${SOURCE_DIR}

#make -j 2 # BUILD EVERYTHING: TODO: BUILD_SUPPORT_SEPARATE.
make TBB  && make

cd BRAINSTools-build/
make test
