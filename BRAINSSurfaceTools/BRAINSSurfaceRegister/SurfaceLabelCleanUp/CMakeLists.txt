StandardBRAINSBuildMacro(NAME SurfaceLabelCleanUp
  ADDITIONAL_SRCS vtkMaskLabel.cxx vtkPolyDataConnectivityIDFilter.cxx
  TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})
