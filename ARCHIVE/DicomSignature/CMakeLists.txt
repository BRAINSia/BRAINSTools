#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(DicomSignature_ITK
  ITKCommon
  ITKIOImageBase
  ITKImageFilterBase
  ITKImageIntensity
  ITKImageStatistics
  ITKStatistics
  ITKTestKernel
  ITKThresholding
)

if(APPLE)
  ## HACK mac does not have openssl headers installed by default.  use "brew install openssl"
  include_directories(/usr/local/opt/openssl/include)
  link_directories(/usr/local/opt/openssl/lib)
endif()



set(ALL_PROGS_LIST
   DicomSignature
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES liblinear ${DicomSignature_ITK_LIBRARIES} BRAINSCommonLib crypto)
endforeach()

### Add the Testing Subdirectory.
#if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
#  add_subdirectory(TestSuite)
#endif()
