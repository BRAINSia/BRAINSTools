#
# Add here command lines to build examples.
#

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_LIST_DIR}/../Testing )

ADD_EXECUTABLE(SmoothingMeshWithScalars_smooth SmoothingMeshWithScalars.cxx)
TARGET_LINK_LIBRARIES(SmoothingMeshWithScalars_smooth ITKQuadEdgeMesh ITKNumerics ITKIO)

ADD_EXECUTABLE(SmoothingMeshWithVectors SmoothingMeshWithVectors.cxx)
TARGET_LINK_LIBRARIES(SmoothingMeshWithVectors ITKQuadEdgeMesh ITKNumerics ITKIO)

