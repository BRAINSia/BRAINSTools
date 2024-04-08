#
# ITK
#

FindITKUtil(ConvertBetweenFileFormats_ITK
  ITKImageIntensity
  ITKTestKernel
  )

if(BUILD_WITH_GUI OR BRAINS3_USE_QT OR BRAINSTools_USE_QT)
  # FINDITKUtil(ITKVtkGlue)
  set(LOCAL_VTKQT GUISupportQt) #vtkGUISupportQt)
  FindVTKUtil(
    InteractionStyle #vtkInteractionStyle
    RenderingCore #vtkRenderingCore
    ImagingSources #vtkImagingSources
    ${LOCAL_VTKQT})
endif()

# This project is intended to be built outside the Insight source tree
project(ConvertBetweenFileFormats)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/castconvertConfigure.h.in
                ${CMAKE_CURRENT_BINARY_DIR}/castconvertConfigure.h @ONLY IMMEDIATE )

set(MODULE_FOLDER "Modules-ConvertBetweenFileFormats")

add_executable(ConvertBetweenFileFormats castconvert.cxx  castconvertDicomScalar.cxx castconvertScalar.cxx  castconvertScalar2D.cxx castconvertScalar3D.cxx castconvertScalar4D.cxx castconvertScalar4DA.cxx castconvertScalar2DA.cxx castconvertScalar3DA.cxx castconvertDicomScalarA.cxx)
set_target_properties(ConvertBetweenFileFormats PROPERTIES FOLDER ${MODULE_FOLDER})

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    )

target_link_libraries(ConvertBetweenFileFormats  ${ConvertBetweenFileFormats_ITK_LIBRARIES} )
if(VTK_FOUND)
  target_link_libraries(ConvertBetweenFileFormats ${VTK_LIBRARIES})
  include_directories( ${ITKApps_SOURCE_DIR}/Auxiliary/vtk )
endif()

install_targets(/bin ConvertBetweenFileFormats)

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
    add_subdirectory(TestSuite)
endif()
