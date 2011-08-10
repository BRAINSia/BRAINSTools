project( GTRACT )
set(LOCAL_PROJECT_NAME GTRACT)

if (BUILD_AGAINST_SLICER3)
  find_package(Slicer3 REQUIRED)
  if(Slicer3_FOUND)
    include(${Slicer3_USE_FILE})
  else(Slicer3_FOUND)
    message(ERROR "Cannot build without a Slicer3 build tree or a Slicer3 installation. Please set Slicer3_DIR. If using a Slicer3 build, Slicer3_DIR needs to be set to the root of the build tree of Slicer3.  If using an installation
   of Slicer3, Slicer3_DIR needs to be set to <SlicerInstallation>/lib/Slicer3.")
  endif(Slicer3_FOUND)

  # Demian W: THESE 3 LINES ARE A HACK SO THE CMakeBuildMacros OF BRAINSTOOLS RECOGNIZE THAT BEING BUILT AS A MODULE IS THE SAME AS BEING BUILT INSIDE OF SLICER
  set(Slicer_SOURCE_DIR ${Slicer3_CMAKE_DIR}/..)
  set(BRAINS_CMAKE_HELPER_DIR ${Slicer3_CMAKE_DIR}/../Applications/CLI/BRAINSTools/BuildScripts)
  set(BRAINSCommonLib_DIR ${Slicer3_DIR}/Applications/CLI/BRAINSTools/BRAINSCommonLib)
endif(BUILD_AGAINST_SLICER3)

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})


option(DEBUGGING_PRINT_IMAGES "Use the DEBUG_IMAGE_VIEWER for debugging" ON)


if(NOT VTK_FOUND)
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
endif(NOT VTK_FOUND)

find_package(BRAINSCommonLib NO_MODULE REQUIRED)
include(${BRAINSCommonLib_USE_FILE})

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

link_libraries(${ITK_LIBRARIES})
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

