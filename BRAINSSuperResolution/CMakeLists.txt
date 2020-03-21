
#
# ITK
#

FindITKUtil(BRAINSSuperResolution_ITK
  ITKAnisotropicSmoothing
  ITKBinaryMathematicalMorphology
  ITKConnectedComponents
  ITKCurvatureFlow
  ITKDistanceMap
  ITKFFT
  ITKImageCompare
  ITKImageFilterBase
  ITKImageFunction
  ITKImageGradient
  ITKImageGrid
  ITKImageIntensity
  ITKImageSources
  ITKImageStatistics
  ITKMathematicalMorphology
  ITKMetricsv4
  ITKOptimizersv4
  ITKRegionGrowing
  ITKRegistrationCommon
  ITKRegistrationMethodsv4
  ITKSmoothing
  ITKStatistics
  ITKTestKernel
  ITKThresholding
  ITKTransform
  RTK
  )


add_subdirectory(GenerateEdgeMap)
add_subdirectory(EdgeBasedWeightedTV)
