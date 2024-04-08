##- project(BRAINSABC)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
find_package(ITK 5.1 REQUIRED)
include(${ITK_USE_FILE})
set(BRAINSABC_ITK_LIBRARIES ${ITK_LIBRARIES})
include_directories( ${ITK_INCLUDE_DIRS})

#-----------------------------------------------------------------------------
# Output directories.
#

###

set(MODULE_FOLDER "Module-BRAINSABC")

add_subdirectory(brainseg)

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  add_subdirectory(TestSuite)
endif()
