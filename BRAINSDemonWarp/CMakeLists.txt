
##- project(BRAINSDemonWarp)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
findITKUtil( BRAINSDemonWarp
  ITKDisplacementField
  ITKIOImageBase
  ITKIOMeta
  ITKImageAdaptors
  ITKImageCompare
  ITKImageCompose
  ITKImageFunction
  ITKImageGrid
  ITKImageIntensity
  ITKImageStatistics
  ITKPDEDeformableRegistration
  ITKRegistrationCommon
  ITKReview
  ITKSmoothing
  ITKSpatialObjects
  ITKTransform
  ITKTestKernel
)

#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSDemonWarpCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSDemonWarpCommonLib.h
  )


set(BRAINSDemonWarpTemplatesLIB_SRCS
         BRAINSDemonWarpTemplates_Vdouble.cxx
         BRAINSDemonWarpTemplates_Vshort.cxx
         BRAINSDemonWarpTemplates_Vushort.cxx
         BRAINSDemonWarpTemplates_int.cxx
         BRAINSDemonWarpTemplates_uint.cxx
         BRAINSDemonWarpTemplates_Vfloat.cxx
         BRAINSDemonWarpTemplates_Vuchar.cxx
         BRAINSDemonWarpTemplates_double.cxx
         BRAINSDemonWarpTemplates_short.cxx
         BRAINSDemonWarpTemplates_ushort.cxx
         BRAINSDemonWarpTemplates_Vint.cxx
         BRAINSDemonWarpTemplates_Vuint.cxx
         BRAINSDemonWarpTemplates_float.cxx
         BRAINSDemonWarpTemplates_uchar.cxx
         BRAINSDemonWarpTemplates.cxx
)

add_library(BRAINSDemonWarpTemplatesLIB STATIC ${BRAINSDemonWarpTemplatesLIB_SRCS})

set(BRAINSDemonWarpTemplatesLIBLibraries BRAINSCommonLib)
DebugImageViewerLibAdditions(BRAINSDemonWarpTemplatesLIBLibraries)

target_link_libraries(BRAINSDemonWarpTemplatesLIB ${BRAINSDemonWarpTemplatesLIBLibraries})

#
# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
#
# XXX When minimum CMake version will be changed to version >=2.8.9, the following code
# could be used instead:
# set_target_properties(BRAINSDemonWarpTemplatesLIB PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set_target_properties(BRAINSDemonWarpTemplatesLIB PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

set(ALL_PROGS_LIST
  BRAINSDemonWarp
  VBRAINSDemonWarp
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSDemonWarpTemplatesLIB )
endforeach()

### Add the Testing Subdirectory.
if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()
