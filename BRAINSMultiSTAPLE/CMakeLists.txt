
FindITKUtil(BRAINSMultiSTAPLE_ITK
  ITKImageCompare
  ITKLabelVoting
  ITKTestKernel
  )
StandardBRAINSBuildMacro( NAME BRAINSMultiSTAPLE
  TARGET_LIBRARIES BRAINSCommonLib ${BRAINSMultiSTAPLE_ITK_LIBRARIES} )

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
    add_subdirectory(TestSuite)
endif()

