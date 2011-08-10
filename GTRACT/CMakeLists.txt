project( GTRACT )
set(LOCAL_PROJECT_NAME GTRACT)

option(DEBUGGING_PRINT_IMAGES "Use the DEBUG_IMAGE_VIEWER for debugging" ON)

#-----------------------------------------------------------------------------
# Output directories.
#

###
# On Visual Studio 8 MS deprecated C. This removes all 1.276E1265 security
# warnings
if(WIN32)
  if(NOT BORLAND)
    if(NOT CYGWIN)
      if(NOT MINGW)
        if(NOT ITK_ENABLE_VISUAL_STUDIO_DEPRECATED_C_WARNINGS)
          add_definitions(
            -D_CRT_FAR_MAPPINGS_NO_DEPRECATE
            -D_CRT_IS_WCTYPE_NO_DEPRECATE
            -D_CRT_MANAGED_FP_NO_DEPRECATE
            -D_CRT_NONSTDC_NO_DEPRECATE
            -D_CRT_SECURE_NO_DEPRECATE
            -D_CRT_SECURE_NO_DEPRECATE_GLOBALS
            -D_CRT_SETERRORMODE_BEEP_SLEEP_NO_DEPRECATE
            -D_CRT_TIME_FUNCTIONS_NO_DEPRECATE
            -D_CRT_VCCLRIT_NO_DEPRECATE
            -D_SCL_SECURE_NO_DEPRECATE
            )
        endif(NOT ITK_ENABLE_VISUAL_STUDIO_DEPRECATED_C_WARNINGS)
      endif(NOT MINGW)
    endif(NOT CYGWIN)
  endif(NOT BORLAND)
endif(WIN32)

if(NOT VTK_FOUND)
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
endif(NOT VTK_FOUND)
link_libraries( vtkGraphics vtkImaging vtkIO vtkFiltering vtkCommon )
#-----------------------------------------------------------------------------
# Configure files with settings for use by the build.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gtractConfigure.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/gtractConfigure.h)

# Definitions for GTRACT Includes
include_directories( ${CMAKE_CURRENT_BINARY_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/Common )

#-----------------------------------------------------------------------------
# Build Generate CLP
#add_subdirectory( Utils )

#-----------------------------------------------------------------------------
# Build GTRACT Tools
add_subdirectory(Common)
add_subdirectory(Cmdline)

