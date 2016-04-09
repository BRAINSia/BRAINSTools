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
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -DANTS_SUPERBUILD:BOOL=OFF
  -DUSE_VTK:BOOL=OFF
  )
if(${PRIMARY_PROJECT_NAME}_USE_QT)
  list(APPEND ${proj}_CMAKE_OPTIONS -DANTS_USE_QT:BOOL=ON)
endif()
### --- End Project specific additions
set(${proj}_REPOSITORY "https://github.com/stnava/ANTs.git")
set(${proj}_GIT_TAG 64864762913c04742643d651985ebc20798c445c) # 20160314
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
