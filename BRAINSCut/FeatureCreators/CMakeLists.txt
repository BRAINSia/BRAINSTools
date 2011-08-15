project(AllFeatureVectors)

##  In many cases sub-projects depending on SlicerExectuion Model
##  that can be built stand alone are combined in larger packages
##  This logic will include SlicerExectionModel only if it
##  has not already been built.

add_subdirectory(DilateImage)
add_subdirectory(DistanceMaps)
add_subdirectory(ErodeImage)
add_subdirectory(FlippedDifference)
add_subdirectory(NeighborhoodMean)
add_subdirectory(NeighborhoodMedian)
add_subdirectory(CannyEdge)
##add_subdirectory(CoOccurrenceMatrix)
add_subdirectory(HammerAttributesCreater)
#add_subdirectory(ImageDistribution)
add_subdirectory(ImageRegionPlotter)
add_subdirectory(JointHistogram)
add_subdirectory(HistogramMatchingFilter)
add_subdirectory(TestImageGeneration)
