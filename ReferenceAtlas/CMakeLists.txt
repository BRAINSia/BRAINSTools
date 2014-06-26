project(ReferenceAtlas)
#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMake)

#-----------------------------------------------------------------------
# Setup locations to find externally maintained test data.
#-----------------------------------------------------------------------
include(ReferenceAtlasExternalData)

# Tell ExternalData commands to transform raw files to content links.
set(ExternalData_LINK_CONTENT MD5)
set(ExternalData_SOURCE_ROOT ${${PROJECT_NAME}_SOURCE_DIR})
include(ExternalData)


#
# In order to test need an XML file

## Configure the standard version
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ExtendedAtlasDefinition.xml.in
  ${ATLAS_INSTALL_DIRECTORY}/ExtendedAtlasDefinition.xml @ONLY IMMEDIATE)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ExtendedAtlasDefinition.xml.in
  DESTINATION ${ATLAS_INSTALL_DIRECTORY}/ USE_SOURCE_PERMISSIONS)

## A macro to wrap getting data downloads
macro(DownloadAndInstallRules GLOBPAT INSTDEST)
  #message(STATUS "YYYYYY ${GLOBPAT}")
  set(AtlasImageFiles "")
  file(GLOB AtlasImageFilesMD5Names ${GLOBPAT})
  #message(STATUS "ABCDE:  ${AtlasImageFilesMD5Names}")
  foreach(CurrentMD5Reference ${AtlasImageFilesMD5Names})
    get_filename_component(ThisImageVarName ${CurrentMD5Reference} NAME)
    string(REPLACE  ".md5" "" CurrentOutFile ${CurrentMD5Reference})
    ExternalData_expand_arguments(${PROJECT_NAME}FetchData ${ThisImageVarName}_VAR DATA{${CurrentOutFile}})
    list(APPEND AtlasImageFiles ${${ThisImageVarName}_VAR})
  endforeach(CurrentMD5Reference ${AtlasImageFilesMD5Names})
  install(FILES ${AtlasImageFiles} DESTINATION ${INSTDEST})
endmacro()

DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/*.nii.gz.md5" ${ATLAS_INSTALL_DIRECTORY})
DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/*.nrrd.md5" ${ATLAS_INSTALL_DIRECTORY})
DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/*.txt.md5" ${ATLAS_INSTALL_DIRECTORY})
DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/*.ctbl.md5" ${ATLAS_INSTALL_DIRECTORY})
DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/*.fcsv.md5" ${ATLAS_INSTALL_DIRECTORY})
DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/*.csv.md5" ${ATLAS_INSTALL_DIRECTORY})

DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/modelFiles/*.gz.md5" ${ATLAS_INSTALL_DIRECTORY}/modelFiles)
DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/probabilityMaps/*.nii.gz.md5" ${ATLAS_INSTALL_DIRECTORY}/probabilityMaps)
DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/spatialImages/*.nii.gz.md5" ${ATLAS_INSTALL_DIRECTORY}/spatialImagesodelFiles)
DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/SurfaceAtlas/*.vtk.md5" ${ATLAS_INSTALL_DIRECTORY}/SurfaceAtlas)
DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/20111119_BCD/*_50Lmks.*.md5" ${ATLAS_INSTALL_DIRECTORY}/20111119_BCD)
DownloadAndInstallRules("${CMAKE_CURRENT_SOURCE_DIR}/T1OnlyModels/*.gz.md5" ${ATLAS_INSTALL_DIRECTORY}/T1OnlyModels)

ExternalData_add_target( ${PROJECT_NAME}FetchData )  # Name of data management target
