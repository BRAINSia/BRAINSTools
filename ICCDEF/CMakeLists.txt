
##- project(iccdefRegistrationNew)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(
  ITKImageCompare
  ITKFFT
  ITKConnectedComponents
  ITKMathematicalMorphology
  ITKBinaryMathematicalMorphology
  ITKRegionGrowing
  ITKPDEDeformableRegistration
  ITKRegionGrowing
  ITKSmoothing
  ITKDeprecated
  ITKRegistrationCommon
  ITKThresholding
  ITKImageFeature
  ITKV3Compatibility
  ITKSmoothing
  ITKLevelSets
  ITKReview
  # ITKDeprecated
  # ITKDisplacementField
  # ITKDistanceMap
  # ITKFFT
  # ITKIOImageBase
  # ITKIOMeta
  # ITKIOSpatialObjects
  # ITKIOTransformBase
  # ITKImageFeature
  # ITKImageFilterBase
  # ITKImageFunction
  # ITKImageGrid
  # ITKImageIntensity
  # ITKImageStatistics
  # ITKPDEDeformableRegistration
  # ITKRegistrationCommon
  # ITKReview
  # ITKSpatialObjects
  # ITKTransform
)

#-------------------------------------------------------------------------------
# OUTPUT DIRECTORIES
#
##- include_directories(${COMMON_BRAINSCOMMONLIB_SOURCE_DIR})
##- link_directories(${COMMON_BRAINSCOMMONLIB_BINARY_DIR})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt/iccdefRegistration_New" CACHE PATH "Install path prefix,
     prepended onto install directories." FORCE)
endif()

set(ALL_PROGS_LIST
  AverageBrainGenerator
  CreateMask
  iccdefRegistration_New
  )

set(ICCDEFLibraries BRAINSCommonLib )
DebugImageViewerLibAdditions(ICCDEFLibraries)

foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES ${ICCDEFLibraries})
endforeach()

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()
