# mkdir -p BRAINSTools-sysclang-RelWithDebInfo && cmake -S BRAINSTools -B BRAINSTools-sysclang-RelWithDebInfo -C BRAINSTools/Utilities/build_configs/Darwin/BRAINSTools-sysclang-RelWithDebInfo.cmake

# set(CMAKE_CXX_STANDARD 14 CACHE STRING "The language standard to use." FORCE)

set(CMAKE_CXX_COMPILER "/Library/Developer/CommandLineTools/usr/bin/c++" CACHE FILEPATH "Initial cache" FORCE)
set(CMAKE_C_COMPILER "/Library/Developer/CommandLineTools/usr/bin/cc" CACHE FILEPATH "Initial cache" FORCE)

set(EXTERNAL_PROJECT_BUILD_TYPE RelWithDebInfo CACHE STRING "The requested build type" FORCE)
set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "The requested build type" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "If examples are built" FORCE)
set(BUILD_TESTING ON CACHE BOOL "If testing is used" FORCE)

#clang 10 new warnings -Wtautological-bitwise-compare
#clang 10 new warnings -Wbitwise-conditional-parentheses
#clang 10 new warnings -Wmisleading-indentation
#clang 10 new warnings -Wc99-designator
#clang 10 new warnings -Wreorder-init-list
#clang 10 new warnings -Wsizeof-array-div
#clang 10 new warnings -Wxor-used-as-pow
#clang 10 new warnings -Wfinal-dtor-non-final-class

set(BRAINSToools_C_OPTIMIZATION_FLAGS
        -glldb

        -fno-omit-frame-pointer
        -fno-optimize-sibling-calls

        -mtune=native -march=native
        -Wtautological-overlap-compare -Wtautological-compare
        -Wrange-loop-analysis
        CACHE STRING "ITK optimation flags for C compiler" FORCE)

set(BRAINSToools_CXX_OPTIMIZATION_FLAGS
        -glldb

        -fno-omit-frame-pointer
        -fno-optimize-sibling-calls

        -mtune=native -march=native
        -Wtautological-overlap-compare -Wtautological-compare
        -Wrange-loop-analysis
        CACHE STRING "ITK optimation flags for CXX compiler" FORCE)

set(BRAINSTools_REQUIRES_VTK OFF CACHE BOOL "bld optional component" FORCE)
set(BRAINSTools_BUILD_DICOM_SUPPORT ON CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSCreateLabelMapFromProbabilityMaps ON CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSMultiModeSegment ON CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSMultiSTAPLE ON CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSMush ON CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSPosteriorToContinuousClass ON CACHE BOOL "bld optional component" FORCE)
set(USE_DWIConvert ON CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSLandmarkInitializer OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSROIAuto OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSResample OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSSnapShotWriter OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSStripRotation OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSTransformConvert OFF CACHE BOOL "bld optional component" FORCE)
set(USE_ConvertBetweenFileFormats OFF CACHE BOOL "bld optional component" FORCE)
set(USE_ImageCalculator OFF CACHE BOOL "bld optional component" FORCE)

set(USE_BRAINSCreateLabelMapFromProbabilityMaps OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSMultiModeSegment OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSMultiSTAPLE OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSMush OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSPosteriorToContinuousClass OFF CACHE BOOL "bld optional component" FORCE)
set(USE_DWIConvert OFF CACHE BOOL "bld optional component" FORCE)

set(USE_GTRACT OFF CACHE BOOL "bld optional component" FORCE)

set(USE_BRAINSTalairach OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSSuperResolution OFF CACHE BOOL "bld optional component" FORCE)
set(USE_DebugImageViewer OFF CACHE BOOL "bld optional component" FORCE)
set(USE_ITKMatlabIO OFF CACHE BOOL "bld optional component" FORCE)
set(USE_BRAINSConstellationDetectorGUI OFF CACHE BOOL "bld optional component" FORCE)


# -- Fixup ITK remote branch info
#set(USE BRAINSTools_ITKv5_GIT_REPOSITORY git@github.com:hjmjohnson/ITK.git to override CACHE STRING "Alternate git repo" FORCE)
#set(USE BRAINSTools_ITKv5_GIT_TAG fix-some-error-pr-request CACHE STRING "Alternate git tag" FORCE)
