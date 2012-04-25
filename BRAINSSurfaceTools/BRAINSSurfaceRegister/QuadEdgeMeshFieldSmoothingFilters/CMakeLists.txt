OPTION(USE_VTK "Use VTK rendering" OFF)
OPTION(GENERATE_PAPER_SCREENSHOTS "Generate VTK screenshots for the report" OFF)
OPTION(GENERATE_REPORT "Generate the PDF report from latex files, source code examples and result screenshots" OFF)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_LIST_DIR}/Source
  ${CMAKE_CURRENT_LIST_DIR}/../QuadEdgeMeshRigidRegistration/Source
  ${CMAKE_CURRENT_LIST_DIR}/../QuadEdgeMeshRigidRegistration/Testing
  )

ADD_SUBDIRECTORY( Source )
#ADD_SUBDIRECTORY( Testing )
ADD_SUBDIRECTORY( Examples )

# IF(GENERATE_REPORT)
#   ADD_SUBDIRECTORY( Documents )
# ENDIF(GENERATE_REPORT)
