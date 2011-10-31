
include_directories(
  ${BRAINSABC_SOURCE_DIR}/brainseg
  ${BRAINSABC_SOURCE_DIR}/common
)
#  ${CMAKE_CURRENT_SOURCE_DIR}/../qhull

set(QHULL_SRCS
  ../qhull/user.c
  ../qhull/global.c
  ../qhull/stat.c
  ../qhull/io.c
  ../qhull/geom2.c
  ../qhull/poly2.c
  ../qhull/merge.c
  ../qhull/qhull.c
  ../qhull/geom.c
  ../qhull/poly.c
  ../qhull/qset.c
  ../qhull/mem.c
  QHullMSTClusteringProcess.cxx
  )

set(BRAINSABCCOMMONLIB_SRCS
  ../common/Log.cxx
  EMSParameters.h
  EMSParameters.cxx
  EMSegmentationFilter.h
  EMSegmentationFilter.txx
  EMSegmentationFilter_float+float.cxx
  AtlasRegistrationMethod_float+float.cxx
  AtlasDefinition.cxx
  filterFloatImages.h
  ExtractSingleLargestRegion.cxx
  ExtractSingleLargestRegion.h
  BRAINSABCUtilities.cxx
  BRAINSABCUtilities.h
)

## Build BRAINSABCCOMMONLIB library
##
add_library(BRAINSABCCOMMONLIB STATIC ${BRAINSABCCOMMONLIB_SRCS})

set(BRAINSABCCOMMONLIBLibraries BRAINSCommonLib ${ITK_LIBRARIES})

DebugImageViewerLibAdditions(BRAINSABCCOMMONLIBLibraries)

target_link_libraries(BRAINSABCCOMMONLIB ${BRAINSABCCOMMONLIBLibraries} )

## Build all the programs
##
set(ALL_PROGS_LIST
  BRAINSABC
  ESLR
  GenerateLabelMapFromProbabilityMap
  )
foreach( prog ${ALL_PROGS_LIST} )
  StandardBRAINSBuildMacro( NAME ${prog} TARGET_LIBRARIES BRAINSABCCOMMONLIB )
endforeach()

if(0)
  StandardBRAINSBuildMacro(NAME BRAINSCleanMask                    ADDITIONAL_SRCS ${BRAINSABCCOMMONLIB_SRCS})

  ## A utility program
  add_executable(StandardizeMaskIntensity StandardizeMaskIntensity.cxx)
  target_link_libraries(StandardizeMaskIntensity ${ITK_LIBRARIES})
endif()
