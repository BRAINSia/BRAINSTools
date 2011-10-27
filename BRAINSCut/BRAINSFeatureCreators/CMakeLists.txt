project(BRAINSFeatureCreaters)
cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)



add_subdirectory(CannyEdge)
add_subdirectory(DilateImage)
add_subdirectory(DistanceMaps)
add_subdirectory(ErodeImage)
add_subdirectory(FlippedDifference)
add_subdirectory(NeighborhoodMean)
add_subdirectory(NeighborhoodMedian)
add_subdirectory(ImageRegionPlotter)
add_subdirectory(JointHistogram)
add_subdirectory(GenerateTestImage)

add_subdirectory(HammerAttributeCreator)

add_subdirectory(CannySegmentationLevelSetImageFilter)
add_subdirectory(CombineLabelToMask)
add_subdirectory(DilateMask)
add_subdirectory(GenerateBrainClippedImage)
add_subdirectory(GenerateCsfClippedFromClassifiedImage)
add_subdirectory(GenerateSummedGradientImage)
add_subdirectory(GradientAnisotropicDiffusionImageFilter)

add_subdirectory(HistogramMatchingFilter)
add_subdirectory(TextureMeasureFilter)

# TODO:  This does not seem to compile? add_subdirectory(CoOccurrenceMatrix)
# TODO:  IS THIS STILL VALID?           add_subdirectory(ImageLabelCombine)
