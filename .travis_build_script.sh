#!/usr/bin/env bash

# -e: fail on error
# -v: show commands
# -x: show expanded commands
set -vex

echo "STARTING BUILD"

PROJ_SOURCE_DIR=${TRAVIS_BUILD_DIR}
PROJ_BUILD_DIR=${HOME}/build
PROJ_INSTALL_DIR=${HOME}/install
BUILD_TYPE=Release

printenv

mkdir -p ${PROJ_BUILD_DIR}
cd ${PROJ_SOURCE_DIR}

cmake \
    -G Ninja \
    -DEXTERNAL_PROJECT_BUILD_TYPE:STRING=${BUILD_TYPE} \
    -DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
    -DCMAKE_CXX_STANDARD:STRING=11 \
    -DUSE_ANTS:BOOL=OFF \
    -DBRAINSTools_REQUIRES_VTK=OFF \
    -DUSE_BRAINSABC:BOOL=OFF \
    -DUSE_BRAINSDWICleanup:BOOL=OFF \
    -DUSE_BRAINSInitializedControlPoints:BOOL=OFF \
    -DUSE_BRAINSLabelStats:BOOL=OFF \
    -DUSE_BRAINSLandmarkInitializer:BOOL=OFF \
    -DUSE_BRAINSROIAuto:BOOL=OFF \
    -DUSE_BRAINSResample:BOOL=OFF \
    -DUSE_BRAINSSnapShotWriter:BOOL=OFF \
    -DUSE_BRAINSStripRotation:BOOL=OFF \
    -DUSE_BRAINSTransformConvert:BOOL=OFF \
    -DUSE_ConvertBetweenFileFormats:BOOL=OFF \
    -DUSE_ImageCalculator:BOOL=OFF \
    -DUSE_ReferenceAtlas:BOOL=OFF \
    -S ${PROJ_SOURCE_DIR} \
    -B ${PROJ_BUILD_DIR}

#make -j 2
cd ${PROJ_BUILD_DIR}
ninja

cd ${PROJ_BUILD_DIR}/BRAINSTools-${BUILD_TYPE}-EP${BUILD_TYPE}-build
ctest -D ExperimentalStart
ctest -D ExperimentalConfigure
ctest -D ExperimentalBuild -j2
ctest -D ExperimentalTest --schedule-random -j2 --output-on-failure
ctest -D ExperimentalSubmit
ninja test
