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
macro(DownloadAndInstallRules UNIQUE_TARGET_NAME GLOBPAT INSTDEST)
  #message(STATUS "YYYYYY ${GLOBPAT}")
  set(AtlasImageFiles "")
  file(GLOB AtlasImageFilesMD5Names ${GLOBPAT})
  #message(STATUS "ABCDE:  ${AtlasImageFilesMD5Names}")
  foreach(CurrentMD5Reference ${AtlasImageFilesMD5Names})
    get_filename_component(ThisImageVarName ${CurrentMD5Reference} NAME)
    string(REPLACE  ".md5" "" CurrentOutFile ${CurrentMD5Reference})
    ExternalData_expand_arguments(${PROJECT_NAME}FetchData ${ThisImageVarName}_VAR DATA{${CurrentOutFile}})
    get_filename_component(${ThisImageVarName}_VAR_REALPATH ${${ThisImageVarName}_VAR} REALPATH)
    get_filename_component(${ThisImageVarName}_VAR_BASENAME ${${ThisImageVarName}_VAR} NAME)

    #message(STATUS "ThisImageVarName_VAR ${ThisImageVarName}_VAR ${${ThisImageVarName}_VAR} ${${ThisImageVarName}_VAR_REALPATH}")
    #message(STATUS "COPY ${${ThisImageVarName}_VAR_REALPATH} ${INSTDEST}/${${ThisImageVarName}_VAR_BASENAME}")

    add_custom_command(
        OUTPUT ${INSTDEST}/${${ThisImageVarName}_VAR_BASENAME}
        COMMAND ${CMAKE_COMMAND} ARGS -E copy ${${ThisImageVarName}_VAR} ${INSTDEST}/${${ThisImageVarName}_VAR_BASENAME}
        DEPENDS ${${ThisImageVarName}_VAR_REALPATH}
        )
    list(APPEND AtlasImageFiles ${INSTDEST}/${${ThisImageVarName}_VAR_BASENAME})

  endforeach(CurrentMD5Reference ${AtlasImageFilesMD5Names})
  #message(STATUS "XXXXX\n${AtlasImageFiles}\nZZZZZZZ")
  add_custom_target(${UNIQUE_TARGET_NAME} DEPENDS ${AtlasImageFiles})

endmacro()

DownloadAndInstallRules(AtlasNIFTI "${CMAKE_CURRENT_SOURCE_DIR}/*.nii.gz.md5" ${ATLAS_INSTALL_DIRECTORY})
DownloadAndInstallRules(AtlasNRRD  "${CMAKE_CURRENT_SOURCE_DIR}/*.nrrd.md5" ${ATLAS_INSTALL_DIRECTORY})
DownloadAndInstallRules(AtlasTxt   "${CMAKE_CURRENT_SOURCE_DIR}/*.txt.md5" ${ATLAS_INSTALL_DIRECTORY})
DownloadAndInstallRules(AtlasCtbl  "${CMAKE_CURRENT_SOURCE_DIR}/*.ctbl.md5" ${ATLAS_INSTALL_DIRECTORY})
DownloadAndInstallRules(AtlasFcsv  "${CMAKE_CURRENT_SOURCE_DIR}/*.fcsv.md5" ${ATLAS_INSTALL_DIRECTORY})
DownloadAndInstallRules(AtlasCsv   "${CMAKE_CURRENT_SOURCE_DIR}/*.csv.md5" ${ATLAS_INSTALL_DIRECTORY})

DownloadAndInstallRules(AtlasModels  "${CMAKE_CURRENT_SOURCE_DIR}/modelFiles/*.gz.md5" ${ATLAS_INSTALL_DIRECTORY}/modelFiles)
DownloadAndInstallRules(AtlasProb    "${CMAKE_CURRENT_SOURCE_DIR}/probabilityMaps/*.nii.gz.md5" ${ATLAS_INSTALL_DIRECTORY}/probabilityMaps)
DownloadAndInstallRules(AtlasSpatial "${CMAKE_CURRENT_SOURCE_DIR}/spatialImages/*.nii.gz.md5" ${ATLAS_INSTALL_DIRECTORY}/spatialImages)
DownloadAndInstallRules(AtlasSurface "${CMAKE_CURRENT_SOURCE_DIR}/SurfaceAtlas/*.vtk.md5" ${ATLAS_INSTALL_DIRECTORY}/SurfaceAtlas)
DownloadAndInstallRules(AtlasBCD     "${CMAKE_CURRENT_SOURCE_DIR}/20111119_BCD/*_50Lmks.*.md5" ${ATLAS_INSTALL_DIRECTORY}/20111119_BCD)
DownloadAndInstallRules(AtlasT1Only  "${CMAKE_CURRENT_SOURCE_DIR}/T1OnlyModels/*.gz.md5" ${ATLAS_INSTALL_DIRECTORY}/T1OnlyModels)

add_custom_target(InstallReferenceAtlas DEPENDS AtlasNIFTI AtlasNRRD AtlasTxt AtlasCtbl AtlasFcsv AtlasCsv AtlasModels AtlasProb AtlasSpatial AtlasSurface AtlasBCD AtlasT1Only)

ExternalData_add_target( ${PROJECT_NAME}FetchData )  # Name of data management target
