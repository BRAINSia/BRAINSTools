
##- project(BRAINSLabelStats)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSLabelStats_ITK ITKImageStatistics ITKIOXML)

#-----------------------------------------------------------------------------
# Output directories.
#

set(ALL_PROGS_LIST
  BRAINSLabelStats
  )

foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib ${BRAINSLabelStats_ITK_LIBRARIES})
endforeach()

#if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
#    add_subdirectory(TestSuite)
#endif()
