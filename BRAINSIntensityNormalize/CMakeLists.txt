##- project(BRAINSIntensityNormalize)

#include_directories(${CMAKE_CURRENT_BINARY_DIR})

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSIntensityNormalize_ITK
        ITKBinaryMathematicalMorphology
        ITKCommon
        ITKConnectedComponents
        ITKIOImageBase
        ITKImageStatistics
        ITKMathematicalMorphology
        ITKOptimizers
        ITKRegionGrowing
        ITKThresholding)

#-----------------------------------------------------------------------------
# Output directories.
#
set(ALL_PROGS_LIST
        BRAINSIntensityNormalize
        )
foreach(prog ${ALL_PROGS_LIST})
    StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

if(0)
    if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
        add_subdirectory(TestSuite)
    endif()
endif()
