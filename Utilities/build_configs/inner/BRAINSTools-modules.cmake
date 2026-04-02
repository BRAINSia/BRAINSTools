# BRAINSTools-modules.cmake
# CMake init-cache file (-C flag) for USE_* module selection.
# Shared by all inner builds regardless of build type.
# Does NOT set CMAKE_BUILD_TYPE — pass that on the cmake command line.
#
# Usage:
#   cmake -C /path/to/BRAINSTools-modules.cmake -DCMAKE_BUILD_TYPE=Release -S ... -B ...

# -----------------------------------------------------------------------
# Build options
# -----------------------------------------------------------------------
set(BUILD_TESTING   ON  CACHE BOOL "Enable CTest" FORCE)
set(BUILD_EXAMPLES  OFF CACHE BOOL "Skip examples" FORCE)

# -----------------------------------------------------------------------
# BRAINSTools feature flags
# -----------------------------------------------------------------------
set(BRAINSTools_REQUIRES_VTK        ON  CACHE BOOL "VTK-dependent tools" FORCE)
set(BRAINSTools_BUILD_DICOM_SUPPORT ON  CACHE BOOL "DICOM I/O (DWIConvert)" FORCE)

# -----------------------------------------------------------------------
# Tools enabled by 3D Slicer
# -----------------------------------------------------------------------
set(USE_BRAINSFit               ON  CACHE BOOL "" FORCE)
set(USE_BRAINSROIAuto           ON  CACHE BOOL "" FORCE)
set(USE_BRAINSResample          ON  CACHE BOOL "" FORCE)
set(USE_BRAINSTransformConvert  ON  CACHE BOOL "" FORCE)
set(USE_DWIConvert              ON  CACHE BOOL "" FORCE)

# -----------------------------------------------------------------------
# Tools required by BRAINSAutoWorkup (BAW)
# -----------------------------------------------------------------------
set(USE_AutoWorkup                              ON  CACHE BOOL "" FORCE)
set(USE_ReferenceAtlas                          ON  CACHE BOOL "" FORCE)
set(USE_ANTS                                    ON  CACHE BOOL "" FORCE)
set(USE_BRAINSABC                               ON  CACHE BOOL "" FORCE)
set(USE_BRAINSConstellationDetector             ON  CACHE BOOL "" FORCE)
set(USE_BRAINSLandmarkInitializer               ON  CACHE BOOL "" FORCE)
set(USE_BRAINSCreateLabelMapFromProbabilityMaps ON  CACHE BOOL "" FORCE)
set(USE_BRAINSSnapShotWriter                    ON  CACHE BOOL "" FORCE)
set(USE_GTRACT                                  ON  CACHE BOOL "" FORCE)

# -----------------------------------------------------------------------
# Additional useful tools
# -----------------------------------------------------------------------
set(USE_ImageCalculator         ON  CACHE BOOL "" FORCE)
set(USE_ConvertBetweenFileFormats ON CACHE BOOL "" FORCE)
set(USE_BRAINSMush              ON  CACHE BOOL "" FORCE)
set(USE_BRAINSMultiModeSegment  ON  CACHE BOOL "" FORCE)
set(USE_BRAINSLabelStats        ON  CACHE BOOL "" FORCE)
set(USE_BRAINSMultiSTAPLE       ON  CACHE BOOL "" FORCE)
set(USE_BRAINSPosteriorToContinuousClass ON CACHE BOOL "" FORCE)
set(USE_BRAINSStripRotation     ON  CACHE BOOL "" FORCE)

# -----------------------------------------------------------------------
# Archived / disabled tools
# -----------------------------------------------------------------------
set(BUILD_ARCHIVE                       ON  CACHE BOOL "Build old tools from ARCHIVE/" FORCE)
set(USE_BRAINSTalairach                ON  CACHE BOOL "" FORCE)
# Qt5 installed via Homebrew at /opt/homebrew/opt/qt@5
set(Qt5_DIR "/opt/homebrew/opt/qt@5/lib/cmake/Qt5" CACHE PATH "Qt5 cmake config directory" FORCE)
set(USE_DebugImageViewer               ON  CACHE BOOL "" FORCE)
set(BRAINS_DEBUG_IMAGE_WRITE           OFF CACHE BOOL "" FORCE)
set(USE_BRAINSInitializedControlPoints ON  CACHE BOOL "" FORCE)  # try enabling
set(USE_BRAINSSuperResolution          ON  CACHE BOOL "" FORCE)  # requires RTK (Module_RTK=ON in ITK)
set(USE_ITKMatlabIO                    OFF CACHE BOOL "" FORCE)  # requires MATLAB SDK
set(USE_BRAINSConstellationDetectorGUI ON  CACHE BOOL "" FORCE)
set(BRAINSTools_USE_QT                 ON  CACHE BOOL "Enable Qt-dependent GUI tools" FORCE)
