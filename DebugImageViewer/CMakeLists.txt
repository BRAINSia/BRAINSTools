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

find_package(Qt5 COMPONENTS Core Widgets Network REQUIRED)

set(DebugImageViewerQT_HDRS
QDebugImageViewerWindow.h
QImageViewerWidget.h
QImageDisplay.h
)

qt5_wrap_cpp(MOCSrc ${DebugImageViewerQT_HDRS} )

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

target_link_libraries(DebugImageViewer PUBLIC
  ${DebugImageViewer_ITK_LIBRARIES}
  VTK::GUISupportQt
  VTK::RenderingOpenGL2
  VTK::InteractionImage
  VTK::InteractionStyle
  VTK::CommonSystem
  Qt5::Core Qt5::Widgets Qt5::Network)

set_target_properties(DebugImageViewer PROPERTIES FOLDER ${MODULE_FOLDER})

###
add_executable(DebugImageViewerClientTest DebugImageViewerClientTest.cxx )

target_link_libraries(DebugImageViewerClientTest PRIVATE
${DebugImageViewer_ITK_LIBRARIES}
${VTK_LIBRARIES}
)

set_target_properties(DebugImageViewerClientTest PROPERTIES FOLDER ${MODULE_FOLDER})

###
add_executable(DebugImageViewerClientTest2 DebugImageViewerClientTest2.cxx )

target_link_libraries(DebugImageViewerClientTest2 PRIVATE
${DebugImageViewer_ITK_LIBRARIES}
${VTK_LIBRARIES}
)

set_target_properties(DebugImageViewerClientTest2 PROPERTIES FOLDER ${MODULE_FOLDER})

###
add_executable(DebugImageViewerVectorClientTest DebugImageViewerVectorClientTest.cxx
)
target_link_libraries(DebugImageViewerVectorClientTest PRIVATE
${DebugImageViewer_ITK_LIBRARIES}
${VTK_LIBRARIES}
)

set_target_properties(DebugImageViewerVectorClientTest PROPERTIES FOLDER ${MODULE_FOLDER})
