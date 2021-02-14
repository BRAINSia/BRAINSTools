
FindITKUtil( DWIConvert_ITK
        ITKIODCMTK
        ITKIORAW
        ITKDCMTK
        ITKNrrdIO
        )

FindVTKUtil(
  IOLegacy #vtkIOLegacy
  IOXML # vtkIOXML
  CommonDataModel # vtkCommonDataModel
  FiltersCore # vtkFiltersCore
  FiltersGeneral # vtkFiltersGeneral
)

option(DEBUGGING_PRINT_IMAGES "Use the DEBUG_IMAGE_VIEWER for debugging" ON)
### A local macro for configuring all these programs consistantly for the different platforms.
###

include_directories(${CMAKE_SOURCE_DIR}/DWIConvert)

set(ALL_PROGS_LIST
  gtractImageConformity               # support DWI interface
  gtractInvertRigidTransform          # No need DWI interface
  gtractTransformToDisplacementField  # No need DWI interface
  gtractInvertBSplineTransform        # No need DWI interface
  gtractResampleCodeImage             # No need DWI interface
  gtractResampleFibers                # No need DWI interface
  compareTractInclusion               # No need DWI interface
  gtractAnisotropyMap                 # No need DWI interface
  gtractTensor                        # No need DWI interface  ???
  gtractCreateGuideFiber              # No need DWI interface
  gtractFiberTracking                 # No need DWI interface
  gtractResampleAnisotropy            # No need DWI interface
  gtractResampleB0                    # support DWI interface
  extractNrrdVectorIndex              # support DWI interface
  gtractAverageBvalues                # support DWI interface
  gtractConcatDwi                     # support DWI interface
  gtractCostFastMarching              # No need DWI interface
  gtractFastMarchingTracking          # No need DWI interface
  gtractCopyImageOrientation          # support DWI interface
  gtractClipAnisotropy                # support DWI interface
  gtractInvertDisplacementField       # No need DWI interface
  gtractResampleDWIInPlace            # support DWI interface
  gtractCoregBvalues                  # No need DWI interface
  gtractCoRegAnatomy                  # support DWI interface
#  gtractCoRegAnatomyRigid
#  gtractCoRegAnatomyBspline
#  gtractFreeTracking
#  gtractStreamlineTracking
#  gtractGuidedTracking
#  gtractGraphSearchTracking
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES GTRACTCommon BRAINSCommonLib ${VTK_LIBRARIES} DWIConvertSupportLib)
endforeach()

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  add_subdirectory(TestSuite)
endif()