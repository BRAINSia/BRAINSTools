set(ALL_PROGS_LIST
   DicomSignature
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES liblinear ${BRAINSContinuousClass_ITK_LIBRARIES} BRAINSCommonLib crypto)
endforeach()

### Add the Testing Subdirectory.
#if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
#  add_subdirectory(TestSuite)
#endif()
