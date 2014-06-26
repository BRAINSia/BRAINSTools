##- project( GTRACTCommonLib )

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
target_link_libraries(GTRACTCommon BRAINSCommonLib ${VTK_LIBRARIES} ${double-conversion_LIBRARIES})

#
# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
#
# XXX When minimum CMake version will be changed to version >=2.8.9, the following code
# could be used instead:
# set_target_properties(GTRACTCommon PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set_target_properties(GTRACTCommon PROPERTIES COMPILE_FLAGS "-fPIC")
endif()
