project(MedicalSlicerViewerGui)

  # Find VTK.
  find_package(VTK REQUIRED)
  if(VTK_FOUND)
    include(${VTK_USE_FILE})
  endif(VTK_FOUND)

  set(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
  set(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
  set(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
  set(DESIRED_QT_VERSION ${VTK_DESIRED_QT_VERSION} CACHE FILEPATH "")

  # Find Qt.
  find_package(Qt4 REQUIRED)
  if(QT_USE_FILE)
    include(${QT_USE_FILE})
  else()
    set(QT_LIBRARIES ${QT_QT_LIBRARY})
  endif()

  set( CMAKE_CURRENT_SOURCE_DIRS
  BRAINSConstellationDetectorGUIPrimary.cxx
  QSliceViewer.cxx
  QLabelList.cxx
  QVTKInteractionCallback.cxx
  QDelLabelDialogs.cxx
  QFileDialogs.cxx
  QHelpDialog.cxx
  )

  set( HEADERS
  include/QSliceViewer.h
  include/QLabelList.h
  include/QVTKInteractionCallback.h
  include/QDelLabelDialogs.h
  include/QFileDialogs.h
  include/QHelpDialog.h
  )

  include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/BRAINSHoughEyeDetector
  ${CMAKE_BINARY_DIR}/src
  ${QT_INCLUDE_DIR}
  ${QT_QTGUI_INCLUDE_DIR}
  ${QT_QTCORE_INCLUDE_DIR}
  )

  QT4_WRAP_CPP( HEADERS_MOC ${HEADERS} )

  add_definitions( ${QT_DEFINITIONS} )

add_library( BRAINSConstellationDetectorGUICOMMONLIB STATIC
  ${CMAKE_CURRENT_SOURCE_DIRS}
  ${HEADERS_MOC}
  )
  target_link_libraries( BRAINSConstellationDetectorGUICOMMONLIB
  QVTK
  ${QT_LIBRARIES}
  BRAINSConstellationDetectorCOMMONLIB
  vtkRendering vtkGraphics vtkIO vtkCommon vtkImaging
  )

## Build all the programs
##
set(ALL_PROGS_LIST
  BRAINSConstellationDetectorGUI
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSConstellationDetectorGUICOMMONLIB )
endforeach()



