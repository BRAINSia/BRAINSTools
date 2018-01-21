
##- project(BRAINSRefacer)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSRefacer_ITK
  ITKTransform
  ITKDistanceMap
  ITKThresholding
  ITKTestKernel
  ITKImageCompare
)

#-----------------------------------------------------------------------------
# Output directories.
#

###
set(CMAKE_CXX_STANDARD 11)

set(ALL_PROGS_LIST
  BRAINSRefacer
  )


foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()


if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
    add_subdirectory(TestSuite)
endif()
