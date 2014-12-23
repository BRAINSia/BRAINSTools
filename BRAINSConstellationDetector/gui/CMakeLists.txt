##- project(MedicalSlicerViewerGui)

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

  if(BUILD_WITH_GUI OR BRAINS3_USE_QT OR BRAINSTools_USE_QT)
    set(LOCAL_VTKQT vtkGUISupportQt)
  endif()
  FindVTKUtil(vtkRenderingFreeType vtkRenderingOpenGL vtkInteractionStyle
    vtkRenderingCore  vtkImagingSources ${LOCAL_VTKQT})

  set( BRAINSConstellationDetectorGUI_SRC
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
  ${BRAINSConstellationDetectorGUI_SRC}
  ${HEADERS_MOC}
  )

#
# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
#
# XXX When minimum CMake version will be changed to version >=2.8.9, the following code
# could be used instead:
# set_target_properties(BRAINSConstellationDetectorGUICOMMONLIB PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set_target_properties(BRAINSConstellationDetectorGUICOMMONLIB PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

if(VTK_USE_QVTK)
  set(LOCAL_VTK_QVTK_LIB QVTK)
endif()
  target_link_libraries( BRAINSConstellationDetectorGUICOMMONLIB
 BRAINSCommonLib
 ${LOCAL_VTK_QVTK_LIB}
 ${QT_LIBRARIES}
 ${BRAINSConstellationDetector_ITK_LIBRARIES}
 # BRAINSConstellationDetectorCOMMONLIB
 ${VTK_LIBRARIES}
  )

## Build all the programs
##
set(ALL_PROGS_LIST
  BRAINSConstellationDetectorGUI
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSConstellationDetectorGUICOMMONLIB )
endforeach()



