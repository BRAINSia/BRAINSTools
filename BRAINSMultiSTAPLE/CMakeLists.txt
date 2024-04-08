
FindITKUtil(BRAINSMultiSTAPLE_ITK
  ITKImageCompare
  ITKLabelVoting
  ITKTestKernel
  )
StandardBRAINSBuildMacro( NAME BRAINSMultiSTAPLE
  TARGET_LIBRARIES BRAINSCommonLib ${BRAINSMultiSTAPLE_ITK_LIBRARIES} )

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
    add_subdirectory(TestSuite)
endif()
