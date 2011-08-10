project(BRAINSABC)
set(LOCAL_PROJECT_NAME BRAINSABC)

#-----------------------------------------------------------------------------
# Output directories.
#
#SETOPTIONALDEBUGIMAGEVIEWER()

###
add_subdirectory(brainseg)

if(1)
  add_subdirectory(TestSuite)
endif(1)

