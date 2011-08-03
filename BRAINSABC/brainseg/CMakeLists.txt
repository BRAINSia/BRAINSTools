
#mark_as_advanced(USE_DEBUG_IMAGE_VIEWER)
#if( USE_DEBUG_IMAGE_VIEWER )
#   if(NOT KWWidgets_SOURCE_DIR)
#     find_package(KWWidgets REQUIRED)
#     include(${KWWidgets_USE_FILE})
#   endif(NOT KWWidgets_SOURCE_DIR)
#   add_definitions(-DUSE_DEBUG_IMAGE_VIEWER)
#   find_path(DEBUG_IMAGE_VIEWER_INCLUDE_DIR DebugImageViewerClient.h ${CMAKE_INSTALL_PREFIX}/include)
#   include_directories(${DEBUG_IMAGE_VIEWER_INCLUDE_DIR})
#endif( USE_DEBUG_IMAGE_VIEWER )

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

set(CLP BRAINSABC)
if(0) # Build against Slicer
  ## Include the Slicer macro for setting up default locations!
  SlicerMacroBuildCLI(
    NAME ${CLP}
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
  )
else()
  SEMMacroBuildCLI(
    NAME ${CLP}
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    ADDITIONAL_SRCS ${BRAINSABCCOMMONLIB_SRCS}
    VERBOSE
    )
endif()


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
