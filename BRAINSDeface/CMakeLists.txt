
##- project(BRAINSDeface)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSDeface_ITK
)

#-----------------------------------------------------------------------------
# Output directories.
#

###
set(ALL_PROGS_LIST
  BRAINSDeface
  )

foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

#if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
#    add_subdirectory(TestSuite)
#endif()
