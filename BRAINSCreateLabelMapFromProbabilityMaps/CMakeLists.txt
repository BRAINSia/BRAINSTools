
#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSCreateLabelMapFromProbabilityMaps_ITK
 ITKImageCompare ITKTestKernel)

StandardBRAINSBuildMacro( NAME BRAINSCreateLabelMapFromProbabilityMaps TARGET_LIBRARIES BRAINSCommonLib ${BRAINSCreateLabelMapFromProbabilityMaps_ITK_LIBRARIES})

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
    add_subdirectory(TestSuite)
endif()
