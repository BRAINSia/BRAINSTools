# cmake -S BRAINSTools -B BRAINSTools-Release -C BRAINSTools-Refactor-Release.cmake
#
# Release superbuild config for Darwin (Apple Silicon / x86_64).
# Key differences from BRAINSTools-Refactor-Debug.cmake:
#   - CMAKE_BUILD_TYPE = Release
#   - No -fno-omit-frame-pointer / -fno-optimize-sibling-calls
#     (those are debug/profiling aids; they slow down all EP builds including FFTW and ITK)
#   - Warning flags still present for compile-time issue detection

set(EXTERNAL_PROJECT_BUILD_TYPE Release CACHE STRING "The requested build type" FORCE)
set(CMAKE_BUILD_TYPE            Release CACHE STRING "The requested build type" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "If examples are built" FORCE)
set(BUILD_TESTING   ON CACHE BOOL "If testing is used" FORCE)

# Optimization flags for EP builds (ITK, FFTW, ANTs, VTK, etc.)
# Omits -fno-omit-frame-pointer and -fno-optimize-sibling-calls intentionally:
# those flags prevent tail-call optimisation and frame-pointer elision, hurting
# performance in tight loops.  They belong only in Debug builds.
set(BRAINSToools_C_OPTIMIZATION_FLAGS
  -mtune=native -march=native
  -Wtautological-overlap-compare -Wtautological-compare
  -Wtautological-bitwise-compare
  -Wbitwise-conditional-parentheses
  -Wrange-loop-analysis
  -Wmisleading-indentation
  -Wc99-designator  -Wreorder-init-list
  -Wsizeof-array-div
  -Wxor-used-as-pow
  -Wfinal-dtor-non-final-class
   CACHE STRING "C compiler flags for EP builds" FORCE)

set(BRAINSToools_CXX_OPTIMIZATION_FLAGS
  -mtune=native -march=native
  -Wtautological-overlap-compare -Wtautological-compare
  -Wtautological-bitwise-compare
  -Wbitwise-conditional-parentheses
  -Wrange-loop-analysis
  -Wmisleading-indentation
  -Wc99-designator  -Wreorder-init-list
  -Wsizeof-array-div
  -Wxor-used-as-pow
  -Wfinal-dtor-non-final-class
   CACHE STRING "C++ compiler flags for EP builds" FORCE)

set(BRAINSTools_REQUIRES_VTK                    ON  CACHE BOOL "bld optional component" FORCE)
set(BRAINSTools_USE_GIT_PROTOCOL                OFF CACHE BOOL "Use https instead of git:// for VTK clone" FORCE)
set(BRAINSTools_BUILD_DICOM_SUPPORT             ON  CACHE BOOL "bld optional component" FORCE)

# Tools enabled by Slicer (mirrors Slicer/SuperBuild.cmake BRAINSTools_slicer_options)
set(USE_BRAINSFit                               ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSROIAuto                           ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSResample                          ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSTransformConvert                  ON  CACHE BOOL "bld optional component" FORCE)
set(USE_DWIConvert                              ON  CACHE BOOL "bld optional component" FORCE) # ON because BRAINSTools_BUILD_DICOM_SUPPORT=ON

# Tools needed by BRAINSAutoWorkup (BAW/workflows, BAW/DWIProcessingWorkflows)
set(USE_AutoWorkup                              ON  CACHE BOOL "bld optional component" FORCE)
set(USE_ReferenceAtlas                          ON  CACHE BOOL "bld optional component" FORCE)
set(USE_ANTS                                    ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSABC                               ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSConstellationDetector             ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSLandmarkInitializer               ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSCreateLabelMapFromProbabilityMaps ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSSnapShotWriter                    ON  CACHE BOOL "bld optional component" FORCE)
set(USE_GTRACT                                  ON  CACHE BOOL "bld optional component" FORCE)

# Other useful tools
set(USE_ImageCalculator                         ON  CACHE BOOL "bld optional component" FORCE)
set(USE_ConvertBetweenFileFormats               ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSMush                              ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSMultiModeSegment                  ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSRefacer                           ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSLabelStats                        ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSMultiSTAPLE                       ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSContinuousClass                   ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSPosteriorToContinuousClass        ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSStripRotation                     ON  CACHE BOOL "bld optional component" FORCE)

# Archived Tools unlikely to be used in the future
set(USE_BRAINSTalairach                         OFF CACHE BOOL "bld optional component" FORCE)
set(USE_DebugImageViewer                        OFF CACHE BOOL "bld optional component" FORCE)
set(BRAINS_DEBUG_IMAGE_WRITE                    OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSInitializedControlPoints          ON  CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSSuperResolution                   ON  CACHE BOOL "bld optional component" FORCE)
set(USE_ITKMatlabIO                             OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSConstellationDetectorGUI          OFF CACHE BOOL "bld optional component" FORCE)

# -- Qt5 (keg-only Homebrew install; used by DebugImageViewer and BRAINSConstellationDetectorGUI)
set(Qt5_DIR "/opt/homebrew/opt/qt@5/lib/cmake/Qt5" CACHE PATH "Qt5 cmake config directory" FORCE)

# -- Fixup ITK remote branch info (uncomment to override)
set(BRAINSTools_ITKv5_GIT_REPOSITORY "https://github.com/hjmjohnson/ITK.git" CACHE STRING "Alternate git repo" FORCE)
set(BRAINSTools_ITKv5_GIT_TAG "fftw-compute-optimized-defaults" CACHE STRING "Alternate git tag" FORCE)
