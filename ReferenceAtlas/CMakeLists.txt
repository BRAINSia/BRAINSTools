project(ReferenceAtlas)
#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_LIST_DIR}/CMake)

#-----------------------------------------------------------------------
# Setup locations to find externally maintained test data.
#-----------------------------------------------------------------------
include(ReferenceAtlasExternalData)

# Tell ExternalData commands to transform raw files to content links.
set(ExternalData_LINK_CONTENT SHA512)
set(ExternalData_SOURCE_ROOT ${${PROJECT_NAME}_SOURCE_DIR})
include(ExternalData)

set(MODULE_FOLDER "Module-ReferenceAtlas")

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
  file(GLOB AtlasImageFilesALGONames ${GLOBPAT})
  #message(STATUS "ABCDE:  ${AtlasImageFilesALGONames}")
  foreach(CurrentALGOReference ${AtlasImageFilesALGONames})
    get_filename_component(ThisImageVarName ${CurrentALGOReference} NAME)
    string(REPLACE  ".md5" "" CurrentOutFile_tmp ${CurrentALGOReference})
    string(REPLACE  ".sha512" "" CurrentOutFile ${CurrentOutFile_tmp})
    ExternalData_expand_arguments(${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} ${ThisImageVarName}_VAR DATA{${CurrentOutFile}})
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

  endforeach()
  #message(STATUS "XXXXX\n${AtlasImageFiles}\nZZZZZZZ")
  #message(STATUS "${UNIQUE_TARGET_NAME} DEPENDS ${AtlasImageFiles}")
  add_custom_target(${UNIQUE_TARGET_NAME} DEPENDS ${AtlasImageFiles})
  set_target_properties(${UNIQUE_TARGET_NAME} PROPERTIES FOLDER ${MODULE_FOLDER})

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
DownloadAndInstallRules(AtlasBCD     "${CMAKE_CURRENT_SOURCE_DIR}/20141004_BCD/*_50Lmks.*.md5" ${ATLAS_INSTALL_DIRECTORY}/20141004_BCD)
DownloadAndInstallRules(AtlasBCDSHA512 "${CMAKE_CURRENT_SOURCE_DIR}/20141004_BCD/*_50Lmks.*.sha512" ${ATLAS_INSTALL_DIRECTORY}/20141004_BCD)
DownloadAndInstallRules(AtlasT1Only  "${CMAKE_CURRENT_SOURCE_DIR}/T1OnlyModels/*.gz.md5" ${ATLAS_INSTALL_DIRECTORY}/T1OnlyModels)

add_custom_target(InstallReferenceAtlas
                  DEPENDS AtlasNIFTI AtlasNRRD AtlasTxt AtlasCtbl AtlasFcsv AtlasCsv AtlasModels AtlasProb AtlasSpatial
                          AtlasSurface AtlasBCD AtlasBCDSHA512 AtlasT1Only)
set_target_properties(InstallReferenceAtlas PROPERTIES FOLDER ${MODULE_FOLDER})
