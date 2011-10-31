project(BRAINSABC)
#-----------------------------------------------------------------------------
# Output directories.
#

###
add_subdirectory(brainseg)

if(BUILD_TESTING)
  add_subdirectory(TestSuite)
endif()

