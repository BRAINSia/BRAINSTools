
##- project(BRAINSMultiModeSegment)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSMultiModeSegment_ITK
  ITKImageIntensity
  ITKBinaryMathematicalMorphology
  ITKConnectedComponents
  ITKMathematicalMorphology
  ITKRegionGrowing
  ITKImageCompare
  ITKThresholding)

#-----------------------------------------------------------------------------
# Output directories.
#

###
set(ALL_PROGS_LIST
  BRAINSMultiModeSegment
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  #  add_subdirectory(TestSuite)
endif()
