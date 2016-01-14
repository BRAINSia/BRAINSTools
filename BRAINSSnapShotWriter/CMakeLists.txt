
##- project(BRAINSSnapShotWriter)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSSnapShotWriter_ITK
  ITKIOImageBase
  ITKImageCompose
  ITKImageFilterBase
  ITKImageFusion
  ITKImageGrid
  ITKImageIntensity
  ITKTestKernel
)

#-----------------------------------------------------------------------------
# Output directories.
#

###
set(ALL_PROGS_LIST
  BRAINSSnapShotWriter
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

#if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
#    add_subdirectory(TestSuite)
#endif()
