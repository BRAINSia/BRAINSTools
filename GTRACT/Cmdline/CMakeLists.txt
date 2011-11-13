option(DEBUGGING_PRINT_IMAGES "Use the DEBUG_IMAGE_VIEWER for debugging" ON)
### A local macro for configuring all these programs consistantly for the different platforms.
###
set(ALL_PROGS_LIST
  gtractImageConformity
  gtractInvertRigidTransform
  gtractTransformToDisplacementField
  gtractInvertBSplineTransform
  gtractResampleCodeImage
  gtractResampleFibers
  compareTractInclusion
  gtractAnisotropyMap
  gtractTensor
  gtractCreateGuideFiber
  gtractFiberTracking
  gtractResampleAnisotropy
  gtractResampleB0
  extractNrrdVectorIndex
  gtractAverageBvalues
  gtractConcatDwi
  gtractCostFastMarching
  gtractFastMarchingTracking
  gtractCopyImageOrientation
  gtractClipAnisotropy
  gtractInvertDisplacementField
  gtractResampleDWIInPlace
  gtractCoregBvalues
  gtractCoRegAnatomy
#  gtractCoRegAnatomyRigid
#  gtractCoRegAnatomyBspline
#  gtractFreeTracking
#  gtractStreamlineTracking
#  gtractGuidedTracking
#  gtractGraphSearchTracking
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES GTRACTCommon )
endforeach()


