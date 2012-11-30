
set(DWIConvert_DEPENDENCIES ITKv4 SlicerExecutionModel)

SlicerMacroCheckExternalProjectDependency(DWIConvert)

set(proj DWIConvert)

set(DWIConvert_REPOSITORY
  ${git_protocol}://github.com/BRAINSia/DWIConvert.git)

#set(DWIConvert_GIT_TAG 3df5666d9ba100e15fa424b0c0af4176f26a0bf6)
set(DWIConvert_GIT_TAG f4a18ade3f84f1be010aa6147119f5004bbd4ff0)

ExternalProject_Add(${proj}
  GIT_REPOSITORY ${DWIConvert_REPOSITORY}
  GIT_TAG ${DWIConvert_GIT_TAG}
  ${slicer_external_update}
  SOURCE_DIR ${proj}
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DBUILD_SHARED_LIBS:BOOL=ON
  -DUSE_SYSTEM_DCMTK:BOOL=ON
  -DITK_DIR=${ITK_DIR}
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  DEPENDS
  ${DWIConvert_DEPENDENCIES}
  )

ExternalProject_Add_Step(${proj} InstallSlicerCMakeLists
  COMMENT "Install simple CMakeList.txt for DWIConvert"
  DEPENDEES download
  DEPENDERS configure
  COMMAND ${CMAKE_COMMAND}
  -E copy ${CMAKE_CURRENT_LIST_DIR}/DWIConvert.CMakeLists.txt
  <SOURCE_DIR>/CMakeLists.txt
  )

set(DWIConvert_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
