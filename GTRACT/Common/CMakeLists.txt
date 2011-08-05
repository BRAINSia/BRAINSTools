project( GTRACTCommonLib )

set(GTRACTCommon_SRC
  algo.cxx
  itkTimeSeriesVersorRigidFilter.cxx
  itkTimeSeriesVersorScaleSkewFilter.cxx
  itkAnatomicalVersorRigidFilter.cxx
  itkAnatomicalBSplineFilter.cxx
  itkInvertBSplineFilter.cxx
  itkTensorToAnisotropyImageFilter.cxx
  itkEigenVectorToColorImageFilter.cxx
  itkComputeDiffusionTensorImageFilter.cxx
  itkGtractImageIO.cxx
  itkGtractParameterIO.cxx
  # itkGtractTransformIO.cxx
  # itkDtiTrackingFilterFibersAndPoints.cxx
)


add_library(GTRACTCommon ${GTRACTCommon_SRC})
target_link_libraries(GTRACTCommon ${ITK_LIBRARIES} ${VTK_LIBRARIES})

if (BUILD_AGAINST_SLICER3)
  slicer3_set_plugins_output_path(GTRACTCommon)
  slicer3_install_plugins(GTRACTCommon)
else(BUILD_AGAINST_SLICER3)
  install(TARGETS
    GTRACTCommon
                   RUNTIME DESTINATION bin
                   LIBRARY DESTINATION lib
                   ARCHIVE DESTINATION lib/static)
endif (BUILD_AGAINST_SLICER3)

