INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_LIST_DIR}/Source
  ${CMAKE_CURRENT_BINARY_DIR}/Source
  ${CMAKE_CURRENT_LIST_DIR}/../QuadEdgeMeshRigidRegistration/Source
  ${CMAKE_CURRENT_LIST_DIR}/../QuadEdgeMeshRigidRegistration/Testing
  )

OPTION(BUILD_REPORTS "Generate articles for the Insight Journal" OFF)

OPTION(USE_TESTING_SET "Run tests with Testing data collection" OFF)
OPTION(USE_TRAINING_SET "Run tests with Training data collection" OFF)

CONFIGURE_FILE(
  ${CMAKE_CURRENT_LIST_DIR}/Source/itkQuadEdgeMeshSphericalDiffeomorphicDemonsRegistrationConfigure.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/Source/itkQuadEdgeMeshSphericalDiffeomorphicDemonsRegistrationConfigure.h
  )

ADD_SUBDIRECTORY( Source )
#ADD_SUBDIRECTORY( Testing )
ADD_SUBDIRECTORY( Examples )

IF( BUILD_REPORTS )
  ADD_SUBDIRECTORY( Documents )
ENDIF( BUILD_REPORTS )
