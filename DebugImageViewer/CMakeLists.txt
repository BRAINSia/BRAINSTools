#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil( DebugImageViewer_ITK
  ITKImageGrid
  ITKImageIntensity
)

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

set(MODULE_FOLDER "Module-DebugImageViewer")

###
add_executable(DebugImageViewer
  ${EXE_SRCS}
)

target_link_libraries(DebugImageViewer
  ${DebugImageViewer_ITK_LIBRARIES}
  ${VTK_LIBRARIES}
  ${QT_LIBRARIES})

set_target_properties(DebugImageViewer PROPERTIES FOLDER ${MODULE_FOLDER})

###
add_executable(DebugImageViewerClientTest DebugImageViewerClientTest.cxx )

target_link_libraries(DebugImageViewerClientTest
${DebugImageViewer_ITK_LIBRARIES}
${VTK_LIBRARIES}
)

set_target_properties(DebugImageViewerClientTest PROPERTIES FOLDER ${MODULE_FOLDER})

###
add_executable(DebugImageViewerClientTest2 DebugImageViewerClientTest2.cxx )

target_link_libraries(DebugImageViewerClientTest2
${DebugImageViewer_ITK_LIBRARIES}
${VTK_LIBRARIES}
)

set_target_properties(DebugImageViewerClientTest2 PROPERTIES FOLDER ${MODULE_FOLDER})

###
add_executable(DebugImageViewerVectorClientTest DebugImageViewerVectorClientTest.cxx
)
target_link_libraries(DebugImageViewerVectorClientTest
${DebugImageViewer_ITK_LIBRARIES}
${VTK_LIBRARIES}
)

set_target_properties(DebugImageViewerVectorClientTest PROPERTIES FOLDER ${MODULE_FOLDER})
