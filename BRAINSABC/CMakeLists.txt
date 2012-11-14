##- project(BRAINSABC)
#-----------------------------------------------------------------------------
# Output directories.
#

###
add_subdirectory(brainseg)

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()

