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
# ITK include dirs are provided transitively by ITK_USE_FILE and target_link_libraries

#-----------------------------------------------------------------------------
# Output directories.
#

###

set(MODULE_FOLDER "Module-BRAINSABC")

add_subdirectory(brainseg)

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  add_subdirectory(TestSuite)
endif()
