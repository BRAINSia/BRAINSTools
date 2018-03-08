


include_directories(
  ${BRAINSTools_SOURCE_DIR}/BRAINSABC/brainseg
  ${BRAINSTools_SOURCE_DIR}/BRAINSABC/common
)

#set(QHULL_SRCS
#  QHullMSTClusteringProcess.cxx
#  )

set(BRAINSABCCOMMONLIB_SRCS
  ../common/Log.cxx
  EMSParameters.h
  EMSParameters.cxx
  EMSegmentationFilter.h
  EMSegmentationFilter.hxx
  EMSegmentationFilter_float+float.cxx
  AtlasRegistrationMethod_float+float.cxx
  AtlasDefinition.cxx
  filterFloatImages.h
  BRAINSABCUtilities.cxx
  BRAINSABCUtilities.h
)

## Build BRAINSABCCOMMONLIB library
##
add_library(BRAINSABCCOMMONLIB STATIC ${BRAINSABCCOMMONLIB_SRCS})
add_dependencies(BRAINSABCCOMMONLIB InstallReferenceAtlas )
set_property(TARGET BRAINSABCCOMMONLIB PROPERTY CXX_STANDARD 11)
set_property(TARGET BRAINSABCCOMMONLIB PROPERTY CXX_STANDARD_REQUIRED ON)
set_target_properties(BRAINSABCCOMMONLIB PROPERTIES FOLDER ${MODULE_FOLDER})

set(BABC_CXX11_Required_Features
      cxx_lambdas
      cxx_strong_enums
      cxx_constexpr
      cxx_auto_type)
target_compile_features(BRAINSABCCOMMONLIB PRIVATE ${BABC_CXX11_Required_Features})


set(BRAINSABCCOMMONLIBLibraries BRAINSCommonLib ${BRAINSABC_ITK_LIBRARIES} ${TBB_IMPORTED_TARGETS})

DebugImageViewerLibAdditions(BRAINSABCCOMMONLIBLibraries)

target_link_libraries(BRAINSABCCOMMONLIB ${BRAINSABCCOMMONLIBLibraries} )

#
# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
#
# XXX When minimum CMake version will be changed to version >=2.8.9, the following code
# could be used instead:
# set_target_properties(BRAINSABCCOMMONLIB PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set_target_properties(BRAINSABCCOMMONLIB PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

## Build all the programs
##
set(ALL_PROGS_LIST
  BRAINSABC
  ESLR
  GenerateLabelMapFromProbabilityMap
  GeneratePurePlugMask
  )
foreach( prog ${ALL_PROGS_LIST} )
  StandardBRAINSBuildMacro( NAME ${prog}
    TARGET_LIBRARIES BRAINSABCCOMMONLIB ${double-conversion_LIBRARIES} )
endforeach()

if(0)
  StandardBRAINSBuildMacro(NAME BRAINSCleanMask                    ADDITIONAL_SRCS ${BRAINSABCCOMMONLIB_SRCS})

  ## A utility program
  add_executable(StandardizeMaskIntensity StandardizeMaskIntensity.cxx)
  target_link_libraries(StandardizeMaskIntensity ${BRAINSABC_ITK_LIBRARIES})
  set_target_properties(StandardizeMaskIntensity PROPERTIES FOLDER ${MODULE_FOLDER})
endif()
