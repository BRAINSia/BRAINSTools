#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
find_package(ITK COMPONENTS
  ITKCommon
  ITKIOImageBase
  ITKImageGrid
  ITKImageIntensity
  ITKVTK
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

set(VTK_LIBRARIES QVTK vtkCommon vtkRendering)

if(NOT QT4_FOUND)
  find_package(Qt4 COMPONENTS QtCore QtGui QtNetwork REQUIRED)
  include(${QT_USE_FILE})
endif()

set(DebugImageViewerQT_HDRS
QDebugImageViewerWindow.h
QImageViewerWidget.h
QImageDisplay.h
)

QT4_WRAP_CPP(MOCSrc ${DebugImageViewerQT_HDRS} )

set(EXE_SRCS
  DebugImageViewer.cxx
  QDebugImageViewerWindow.cxx
  QImageDisplay.cxx
  QImageViewerWidget.cxx
  vtkKWImage.cxx
  ${MOCSrc}
)

add_executable(DebugImageViewer
  ${EXE_SRCS}
)

target_link_libraries(DebugImageViewer
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
  ${QT_LIBRARIES})

add_executable(DebugImageViewerClientTest DebugImageViewerClientTest.cxx )

target_link_libraries(DebugImageViewerClientTest
${ITK_LIBRARIES}
${VTK_LIBRARIES}
)

add_executable(DebugImageViewerClientTest2 DebugImageViewerClientTest2.cxx )

target_link_libraries(DebugImageViewerClientTest2
${ITK_LIBRARIES}
${VTK_LIBRARIES}
)

add_executable(DebugImageViewerVectorClientTest DebugImageViewerVectorClientTest.cxx
)
target_link_libraries(DebugImageViewerVectorClientTest
${ITK_LIBRARIES}
${VTK_LIBRARIES}
)
