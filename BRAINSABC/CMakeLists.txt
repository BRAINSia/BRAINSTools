project(BRAINSABC)
set(LOCAL_PROJECT_NAME BRAINSABC)

#-----------------------------------------------------------------------------
# Output directories.
#
#SETOPTIONALDEBUGIMAGEVIEWER()

###
add_subdirectory(brainseg)

if(BUILD_TESTING)
  add_subdirectory(TestSuite)
endif()

