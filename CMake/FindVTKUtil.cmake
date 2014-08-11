#
# find and incorporate the VTK library
macro(FindVTKUtil)
  if(BRAINSTools_REQUIRES_VTK)
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
  endif()
endmacro()