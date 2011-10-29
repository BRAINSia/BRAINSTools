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

add_library(GTRACTCommon STATIC ${GTRACTCommon_SRC})
target_link_libraries(GTRACTCommon BRAINSCommonLib ${VTK_LIBRARIES})

