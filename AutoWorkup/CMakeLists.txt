##- project(AutoWorkup)

#-----------------------------------------------------------------------------
# Output directories.
#
if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
    add_subdirectory(TestSuite)
endif()
