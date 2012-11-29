# This project is intended to be built outside the Insight source tree
project(ConvertBetweenFileFormats)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/castconvertConfigure.h.in
                ${CMAKE_CURRENT_BINARY_DIR}/castconvertConfigure.h @ONLY IMMEDIATE )

add_executable(ConvertBetweenFileFormats castconvert.cxx  castconvertDicomScalar.cxx castconvertScalar.cxx  castconvertScalar2D.cxx castconvertScalar3D.cxx castconvertScalar4D.cxx castconvertScalar4DA.cxx castconvertScalar2DA.cxx castconvertScalar3DA.cxx castconvertDicomScalarA.cxx)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    )

target_link_libraries(ConvertBetweenFileFormats  ${ITK_LIBRARIES} )
if(VTK_FOUND)
  target_link_libraries(ConvertBetweenFileFormats vtkIO vtkCommon vtkFiltering vtkImaging )
  include_directories( ${ITKApps_SOURCE_DIR}/Auxiliary/vtk )
endif(VTK_FOUND)

install_targets(/bin ConvertBetweenFileFormats)

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
    add_subdirectory(TestSuite)
endif()

