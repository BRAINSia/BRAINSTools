
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

StandardBRAINSBuildMacro(NAME BRAINSABC ADDITIONAL_SRCS ${BRAINSABCCOMMONLIB_SRCS})


if(0)

set(${CLP}_SOURCE BRAINSABC.cxx )
set(BRAINSABC_PRIMARY_SOURCE ${BRAINSABC_SOURCE_DIR}/brainseg/BRAINSABCPrimary.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} BRAINSABCPrimary.xml "${${CLP}_SOURCE}" "${BRAINSABC_PRIMARY_SOURCE}" "BRAINSABCPrimary" "")
add_dependencies(BRAINSABCCOMMONLIB BRAINSABCPrimaryCLP.h) ## Needed to create proper library dependancies when building parallel
target_link_libraries( ${CLP} BRAINSABCCOMMONLIB)




set(CLP ESLR)
set(${CLP}_SOURCE ${CLP}.cxx )
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} ${CLP}.xml "${${CLP}_SOURCE}" "" "main" "")
target_link_libraries( ${CLP} BRAINSABCCOMMONLIB )

set(CLP GenerateLabelMapFromProbabilityMap)
set(${CLP}_SOURCE ${CLP}.cxx )
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} ${CLP}.xml "${${CLP}_SOURCE}" "" "main" "")
target_link_libraries( ${CLP} ${ITK_LIBRARIES})

## A utility program
add_executable(StandardizeMaskIntensity StandardizeMaskIntensity.cxx)
target_link_libraries(StandardizeMaskIntensity ${ITK_LIBRARIES})

set(CLP BRAINSCleanMask)
set(${CLP}_SOURCE ${CLP}.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} ${CLP}.xml "${${CLP}_SOURCE}" "" "main" "")
target_link_libraries( ${CLP} ${ITK_LIBRARIES})

endif()
