#
# find and incorporate the VTK library
macro(FindVTKUtil)
  if(NOT BRAINSTools_REQUIRES_VTK)
    message( FATAL_ERROR "You have requested the FindVTKUtil macro, but "
                         "the requesting module is not listed as requiring vtk "
                         "under Common.cmake as \"set(BRAINSTools_REQUIRES_VTK ON)\" "
                         "please add the requesting module to the list in Common.cmake" )
  endif()

  #  message("VTK_DIR:${VTK_DIR}")
  find_package(VTK COMPONENTS
    vtkCommonSystem
    vtkCommonCore
    vtkCommonSystem
    vtkCommonMath
    vtkCommonMisc
    vtkCommonTransforms
    ${ARGN}
    REQUIRED)
  include(${VTK_USE_FILE})
  #  message("VTK_USE_FILE:${VTK_USE_FILE}")
  #  message("VTK_INCLUDE_DIRS:${VTK_INCLUDE_DIRS}")
  include_directories(${VTK_INCLUDE_DIRS})
  link_directories(${VTK_LIBRARY_DIRS})
endmacro()
