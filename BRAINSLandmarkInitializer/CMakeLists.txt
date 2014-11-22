
##- project(BRAINSLandmarkInitializer)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSLandmarkInitializer_ITK
  ITKCommon
  ITKDistanceMap
  ITKImageCompare
  ITKRegistrationCommon
  ITKThresholding
  ITKTransform
)

#-----------------------------------------------------------------------------
# Output directories.
#
#
set(ALL_PROGS_LIST
  BRAINSLandmarkInitializer
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

#if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
#    add_subdirectory(TestSuite)
#endif()
