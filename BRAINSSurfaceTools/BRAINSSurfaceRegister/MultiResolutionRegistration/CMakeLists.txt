CONFIGURE_FILE(
  ${CMAKE_CURRENT_LIST_DIR}/../QuadEdgeMeshSphericalDiffeomorphicDemonsRegistration/Source/itkQuadEdgeMeshSphericalDiffeomorphicDemonsRegistrationConfigure.h.in
  ${CMAKE_CURRENT_LIST_DIR}/itkQuadEdgeMeshSphericalDiffeomorphicDemonsRegistrationConfigure.h
  )

StandardBRAINSBuildMacro(NAME MultiResolutionRegistration TARGET_LIBRARIES ${ITK_LIBRARIES})
