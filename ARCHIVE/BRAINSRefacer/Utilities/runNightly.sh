#!/bin/bash
# Script to run the Test Suite for BRAINSRefacer nightly
# Change the paths of the directories below
# This needs to be copied to and run from a location outside of the source directory.

BASE_DIR=~/code/brains/alexProject
SOURCE_DIR=${BASE_DIR}/BRAINSTools
BUILD_DIR=${BASE_DIR}/nightlyBuild
EXTERNAL_SOURCES_DIR=~/code/brains/nightly/build


if [[ ! -d ${BASE_DIR} ]]; then
  echo "base directory doesn't exist";
  exit -1;
fi

echo "updating source repo"
if [[ ! -d ${SOURCE_DIR} ]]; then
  echo "source directory doesn't exist. Cloning.";
  git clone https://github.com/BRAINSia/BRAINSTools.git;
  cd ${SOURCE_DIR}
  git checkout BRAINSRefacer
else
  cd ${SOURCE_DIR};
  git fetch origin;
  git checkout master;
  git branch -D BRAINSRefacer
  git checkout BRAINSRefacer
fi

echo "changing to build directory"
if [[ ! -d ${BUILD_DIR} ]]; then
  echo "build directory doesnt exist. Creating ..."
  mkdir -p ${BUILD_DIR};
fi

cd ${BUILD_DIR};

echo "Configuring with cmake"
/usr/local/bin/cmake ${SOURCE_DIR} -DBRAINSTools_SUPERBUILD:BOOL=OFF -DSlicerExecutionModel_DIR:PATH=${EXTERNAL_SOURCES_DIR}/SlicerExecutionModel-build -DCMAKE_MODULE_PATH:PATH=${SOURCE_DIR}/CMake -DUSE_BRAINSRefacer=ON -DUSE_ANTS=OFF -DUSE_AutoWorkup=OFF -DUSE_BRAINSABC=OFF -DUSE_BRAINSConstellationDetector=OFF -DUSE_BRAINSDWICleanup=OFF -DUSE_BRAINSFit=OFF -DUSE_BRAINSInitializedControlPoints=OFF -DUSE_BRAINSLabelStats=OFF -DUSE_BRAINSLandmarkInitializer=OFF -DUSE_BRAINSROIAuto=OFF -DUSE_BRAINSResample=OFF -DUSE_BRAINSSnapShotWriter=OFF -DUSE_BRAINSStripRotation=OFF -DUSE_BRAINSTransformConvert=OFF -DUSE_DWIConvert=OFF -DUSE_ConvertBetweenFileFormats=OFF -DUSE_ImageCalculator=OFF -DUSE_ReferenceAtlas=OFF

echo "Building";
make -j24 -k;
make -j24 -k;

echo "Testing"
/usr/local/bin/ctest -j12 -D NightlyStart;
/usr/local/bin/ctest -j12 -D NightlyConfigure;
/usr/local/bin/ctest -j12 -D NightlyBuild;
/usr/local/bin/ctest -j12 -D NightlyTest;
/usr/local/bin/ctest -j12 -D NightlySubmit;
