##- project(MedicalSlicerViewerGui)

  find_package(Qt5 COMPONENTS Core Widgets REQUIRED)

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

  qt5_wrap_cpp( HEADERS_MOC ${HEADERS} )

add_library( BRAINSConstellationDetectorGUICOMMONLIB STATIC
  ${BRAINSConstellationDetectorGUI_SRC}
  ${HEADERS_MOC}
  )
set_target_properties(BRAINSConstellationDetectorGUICOMMONLIB PROPERTIES FOLDER ${MODULE_FOLDER})

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

  target_link_libraries( BRAINSConstellationDetectorGUICOMMONLIB
 BRAINSCommonLib
 VTK::GUISupportQt
 Qt5::Core Qt5::Widgets
 ${BRAINSConstellationDetector_ITK_LIBRARIES}
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
