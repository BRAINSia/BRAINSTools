set(proj        ANTs) #This local name

set(${proj}_DEPENDENCIES ITKv4 SlicerExecutionModel )

if(${PROJECT_NAME}_BUILD_DICOM_SUPPORT)
  list(APPEND ${proj}_DEPENDENCIES DCMTK)
endif()

# Set dependency list
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

### --- Project specific additions here
set(${proj}_CMAKE_OPTIONS
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
  -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
  -DUSE_SYSTEM_ITK:BOOL=ON
  -DUSE_SYSTEM_SlicerExecutionModel:BOOL=ON
  -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
  -DITK_DIR:PATH=${ITK_DIR}
  -DUSE_VTK:BOOL=OFF
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -DANTS_SUPERBUILD:BOOL=OFF
  -DBUILD_ALL_ANTS_APPS:BOOL=ON #Perhaps turn this to OFF
  #  -DANTS_BUILD_DenoiseImage=ON
  #  -DANTS_BUILD_antsRegistration=ON
  #  -DANTS_BUILD_antsJointFusion=ON
  )
if(${PRIMARY_PROJECT_NAME}_USE_QT)
  list(APPEND ${proj}_CMAKE_OPTIONS -DANTS_USE_QT:BOOL=ON)
endif()
### --- End Project specific additions
set(${proj}_REPOSITORY "https://github.com/stnava/ANTs.git")
set(${proj}_GIT_TAG bd9a952048a23c21f3a28d8ba2b9e81a96a4496f) # 20160824
ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  GIT_REPOSITORY ${${proj}_REPOSITORY}
  GIT_TAG ${${proj}_GIT_TAG}
  SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
  BINARY_DIR ${proj}-build
  LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
  LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
  LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
  LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
  ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS -Wno-dev --no-warn-unused-cli
  CMAKE_CACHE_ARGS
    ${${proj}_CMAKE_OPTIONS}
  INSTALL_COMMAND ""
  DEPENDS
  ${${proj}_DEPENDENCIES}
  )

set(${proj}_SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj})
set(${proj}_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib)

mark_as_superbuild(
  VARS ${proj}_SOURCE_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
