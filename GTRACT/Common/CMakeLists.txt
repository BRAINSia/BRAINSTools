##- project( GTRACTCommonLib )

set(GTRACTCommon_SRC
  algo.cxx
  itkAnatomicalVersorRigidFilter.cxx
  itkAnatomicalBSplineFilter.cxx
  itkInvertBSplineFilter.cxx
  itkTensorToAnisotropyImageFilter.cxx
  itkEigenVectorToColorImageFilter.cxx
  itkComputeDiffusionTensorImageFilter.cxx
  itkGtractImageIO.cxx
  itkGtractParameterIO.cxx
)

add_library(GTRACTCommon STATIC ${GTRACTCommon_SRC})
target_include_directories(GTRACTCommon PUBLIC
  ${GTRACT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(GTRACTCommon PUBLIC BRAINSCommonLib ${VTK_LIBRARIES} ${double-conversion_LIBRARIES})

# Suppress MSVC deprecation warnings for all GTRACT targets (propagated via PUBLIC)
if(WIN32 AND NOT BORLAND AND NOT CYGWIN AND NOT MINGW
   AND NOT ITK_ENABLE_VISUAL_STUDIO_DEPRECATED_C_WARNINGS)
  target_compile_definitions(GTRACTCommon PUBLIC
    _CRT_FAR_MAPPINGS_NO_DEPRECATE
    _CRT_IS_WCTYPE_NO_DEPRECATE
    _CRT_MANAGED_FP_NO_DEPRECATE
    _CRT_NONSTDC_NO_DEPRECATE
    _CRT_SECURE_NO_DEPRECATE
    _CRT_SECURE_NO_DEPRECATE_GLOBALS
    _CRT_SETERRORMODE_BEEP_SLEEP_NO_DEPRECATE
    _CRT_TIME_FUNCTIONS_NO_DEPRECATE
    _CRT_VCCLRIT_NO_DEPRECATE
    _SCL_SECURE_NO_DEPRECATE
  )
endif()

set_target_properties(GTRACTCommon PROPERTIES FOLDER ${MODULE_FOLDER})

if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set_target_properties(GTRACTCommon PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
