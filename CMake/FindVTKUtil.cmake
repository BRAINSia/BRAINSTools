#
# find and incorporate the VTK library
macro(FindVTKUtil)
  set(VTK_VERSION_MIN 9.0.1)
  if(NOT ${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK)
    message( FATAL_ERROR "You have requested the FindVTKUtil macro, but "
                         "the requesting module is not listed as requiring vtk "
                         "under Common.cmake as \"set(${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK ON)\" "
                         "please add the requesting module to the list in Common.cmake" )
  endif()

  #  message("VTK_DIR:${VTK_DIR}")
  set(VTK_COMMON_COMPONENTS
    CommonSystem #vtkCommonSystem
    CommonCore #vtkCommonCore
    CommonSystem #vtkCommonSystem
    CommonMath #vtkCommonMath
    CommonMisc #vtkCommonMisc
    CommonTransforms #vtkCommonTransforms
    ${ARGN}
  )
  find_package(VTK ${VTK_VERSION_MIN} COMPONENTS ${VTK_COMMON_COMPONENTS} REQUIRED)

  ## Paradigm specified by http://www.vtk.org/Wiki/VTK/Build_System_Migration#How_Implementation_Modules_Are_Initialized
  include_directories(${VTK_INCLUDE_DIRS})

#  if(TARGET vtkRenderingVolumeOpenGL)
#    message(STATUS "Building optional volume rendering component")
#    find_package(VTK ${VTK_VERSION_MIN} COMPONENTS ${VTK_COMMON_COMPONENTS} vtkRenderingVolumeOpenGL)
##  endif()
#  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${VTK_DEFINITIONS})

  ## Paradigm specified by http://www.vtk.org/Wiki/VTK/Build_System_Migration#How_Implementation_Modules_Are_Initialized
#  include_directories(${VTK_INCLUDE_DIRS})
#  if(TARGET vtkRenderingVolumeOpenGL2)
#    message(STATUS "Building optional volume rendering component")
#    find_package(VTK ${VTK_VERSION_MIN} COMPONENTS ${VTK_COMMON_COMPONENTS} vtkRenderingVolumeOpenGL2)
#  endif()
#  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${VTK_DEFINITIONS})

  find_package(VTK ${VTK_VERSION_MIN} REQUIRED) ## HACK: VTK should be minimized.  This is maximizing it's use to all modules.

  # No longer needed; warns or errors depending on the version requested when
  # finding VTK.
  #include(${VTK_USE_FILE})

  #  message("VTK_USE_FILE:${VTK_USE_FILE}")
  #  message("VTK_INCLUDE_DIRS:${VTK_INCLUDE_DIRS}")
  include_directories(${VTK_INCLUDE_DIRS})
  link_directories(${VTK_LIBRARY_DIRS})
endmacro()
