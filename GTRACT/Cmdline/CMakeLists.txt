option(DEBUGGING_PRINT_IMAGES "Use the DEBUG_IMAGE_VIEWER for debugging" ON)
### A local macro for configuring all these programs consistantly for the different platforms.
macro(ADDGTRACTPROG PROGNAME)
###
if(0) # Build against Slicer
  ## Include the Slicer macro for setting up default locations!
  SlicerMacroBuildCLI(
    NAME ${PROGNAME}
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES GTRACTCommon ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
  )
else()
  SEMMacroBuildCLI(
    NAME ${PROGNAME}
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES GTRACTCommon ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
    )
endif()
endmacro(ADDGTRACTPROG PROGNAME)

ADDGTRACTPROG(gtractImageConformity)
ADDGTRACTPROG(gtractInvertRigidTransform)
ADDGTRACTPROG(gtractTransformToDeformationField)
ADDGTRACTPROG(gtractInvertBSplineTransform)
ADDGTRACTPROG(gtractResampleCodeImage)
ADDGTRACTPROG(gtractResampleFibers)
ADDGTRACTPROG(compareTractInclusion)
ADDGTRACTPROG(gtractAnisotropyMap)
ADDGTRACTPROG(gtractTensor)
ADDGTRACTPROG(gtractCreateGuideFiber)
ADDGTRACTPROG(gtractFiberTracking)
ADDGTRACTPROG(gtractResampleAnisotropy)
ADDGTRACTPROG(gtractResampleB0)
ADDGTRACTPROG(extractNrrdVectorIndex)
ADDGTRACTPROG(gtractAverageBvalues)
ADDGTRACTPROG(gtractConcatDwi)
ADDGTRACTPROG(gtractCostFastMarching)
ADDGTRACTPROG(gtractFastMarchingTracking)
ADDGTRACTPROG(gtractCopyImageOrientation)
ADDGTRACTPROG(gtractClipAnisotropy)
ADDGTRACTPROG(gtractInvertDeformationField)
ADDGTRACTPROG(gtractResampleDWIInPlace)
ADDGTRACTPROG(gtractCoregBvalues)
ADDGTRACTPROG(gtractCoRegAnatomy)

#ADDGTRACTPROG(gtractCoRegAnatomyRigid)
#ADDGTRACTPROG(gtractCoRegAnatomyBspline)
#ADDGTRACTPROG(gtractFreeTracking)
#ADDGTRACTPROG(gtractStreamlineTracking)
#ADDGTRACTPROG(gtractGuidedTracking)
#ADDGTRACTPROG(gtractGraphSearchTracking)

