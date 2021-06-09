set(proj        ANTs) #This local name

set(${proj}_DEPENDENCIES ITKv5 SlicerExecutionModel )

if(${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT)
  list(APPEND ${proj}_DEPENDENCIES DCMTK)
endif()
if(${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK)
  list(APPEND ${proj}_DEPENDENCIES VTK)
endif()


# Set dependency list
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)



### --- Project specific additions here
set(${proj}_CMAKE_OPTIONS
  -DUSE_SYSTEM_ITK:BOOL=ON
  -DITK_DIR:PATH=${ITK_DIR}
  -DUSE_VTK:BOOL=OFF
  -DUSE_SYSTEM_SlicerExecutionModel:BOOL=ON
  -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}

  -DANTS_SUPERBUILD:BOOL=OFF
  -DBUILD_ALL_ANTS_APPS:BOOL=OFF #Perhaps turn this to OFF
  -DANTS_BUILD_WITH_CCACHE:BOOL=OFF # This does not work well
  -DANTS_BUILD_antsAffineInitializer:BOOL=OFF
  -DANTS_BUILD_antsJointFusion:BOOL=ON
  -DANTS_BUILD_DenoiseImage:BOOL=ON
  -DANTS_BUILD_SurfaceBasedSmoothing:BOOL=OFF
  -DANTS_BUILD_ThresholdImage:BOOL=OFF
  -DANTS_BUILD_ResampleImage:BOOL=OFF
  -DANTS_BUILD_sccan:BOOL=OFF
  -DANTS_BUILD_N4BiasFieldCorrection:BOOL=OFF
  -DANTS_BUILD_N3BiasFieldCorrection:BOOL=OFF
  -DANTS_BUILD_KellyKapowski:BOOL=OFF
  -DANTS_BUILD_antsRegistration:BOOL=ON
  -DANTS_BUILD_antsMotionCorrStats:BOOL=OFF
  -DANTS_BUILD_antsMotionCorr:BOOL=OFF
  -DANTS_BUILD_antsApplyTransforms:BOOL=OFF
  -DANTS_BUILD_LabelGeometryMeasures:BOOL=OFF
  -DANTS_BUILD_LabelClustersUniquely:BOOL=OFF
  -DANTS_BUILD_Atropos:BOOL=OFF
  -DANTS_BUILD_antsApplyTransformsToPoints:BOOL=OFF

  -DANTS_BUILD_antsAI:BOOL=OFF
  -DANTS_BUILD_antsJointTensorFusion:BOOL=OFF
  -DANTS_BUILD_ImageMath:BOOL=OFF
  -DANTS_BUILD_iMath:BOOL=OFF
  -DANTS_BUILD_ANTS:BOOL=OFF
  -DANTS_BUILD_ANTSJacobian:BOOL=OFF
  -DANTS_BUILD_CreateJacobianDeterminantImage:BOOL=OFF
  -DANTS_BUILD_PrintHeader:BOOL=OFF
  -DANTS_BUILD_ResetDirection:BOOL=OFF
  -DANTS_BUILD_ANTSUseLandmarkImagesToGetAffineTransform:BOOL=OFF
  -DANTS_BUILD_ANTSUseLandmarkImagesToGetBSplineDisplacementField:BOOL=OFF
  -DANTS_BUILD_ANTSUseDeformationFieldToGetAffineTransform:BOOL=OFF
  -DANTS_BUILD_antsLandmarkBasedTransformInitializer:BOOL=OFF
  -DANTS_BUILD_LaplacianThickness:BOOL=OFF
  -DANTS_BUILD_SetOrigin:BOOL=OFF
  -DANTS_BUILD_SetSpacing:BOOL=OFF
  -DANTS_BUILD_SetDirectionByMatrix:BOOL=OFF
  -DANTS_BUILD_SurfaceCurvature:BOOL=OFF
  -DANTS_BUILD_ConvertScalarImageToRGB:BOOL=OFF
  -DANTS_BUILD_CreateWarpedGridImage:BOOL=OFF
  -DANTS_BUILD_MeasureImageSimilarity:BOOL=OFF
  -DANTS_BUILD_ConvertToJpg:BOOL=OFF
  -DANTS_BUILD_ConvertImage:BOOL=OFF
  -DANTS_BUILD_ConvertImagePixelType:BOOL=OFF
  -DANTS_BUILD_ConvertInputImagePixelTypeToFloat:BOOL=OFF
  -DANTS_BUILD_FitBSplineToPoints:BOOL=OFF
  -DANTS_BUILD_AverageTensorImages:BOOL=OFF
  -DANTS_BUILD_ImageSetStatistics:BOOL=OFF
  -DANTS_BUILD_MultiplyImages:BOOL=OFF
  -DANTS_BUILD_SimulateDisplacementField:BOOL=OFF
  -DANTS_BUILD_SmoothDisplacementField:BOOL=OFF
  -DANTS_BUILD_SmoothImage:BOOL=OFF
  -DANTS_BUILD_ClusterImageStatistics:BOOL=OFF
  -DANTS_BUILD_LabelOverlapMeasures:BOOL=OFF
  -DANTS_BUILD_LesionFilling:BOOL=OFF
  -DANTS_BUILD_MeasureMinMaxMean:BOOL=OFF
  -DANTS_BUILD_WarpImageMultiTransform:BOOL=OFF
  -DANTS_BUILD_ComposeMultiTransform:BOOL=OFF
  -DANTS_BUILD_MemoryTest:BOOL=OFF
  -DANTS_BUILD_PermuteFlipImageOrientationAxes:BOOL=OFF
  -DANTS_BUILD_ImageCompare:BOOL=OFF
  -DANTS_BUILD_ResampleImageBySpacing:BOOL=OFF
  -DANTS_BUILD_CopyImageHeaderInformation:BOOL=OFF
  -DANTS_BUILD_WarpTimeSeriesImageMultiTransform:BOOL=OFF
  -DANTS_BUILD_ExtractSliceFromImage:BOOL=OFF
  -DANTS_BUILD_ExtractRegionFromImage:BOOL=OFF
  -DANTS_BUILD_ExtractRegionFromImageByMask:BOOL=OFF
  -DANTS_BUILD_PasteImageIntoImage:BOOL=OFF
  -DANTS_BUILD_TileImages:BOOL=OFF
  -DANTS_BUILD_CreateTiledMosaic:BOOL=OFF
  -DANTS_BUILD_CreateImage:BOOL=OFF
  -DANTS_BUILD_NonLocalSuperResolution:BOOL=OFF
  -DANTS_BUILD_WarpTensorImageMultiTransform:BOOL=OFF
  -DANTS_BUILD_ReorientTensorImage:BOOL=OFF
  -DANTS_BUILD_RebaseTensorImage:BOOL=OFF
  -DANTS_BUILD_KellySlater:BOOL=OFF
  -DANTS_BUILD_CreateDTICohort:BOOL=OFF
  -DANTS_BUILD_antsAlignOrigin:BOOL=OFF
  -DANTS_BUILD_antsMotionCorrDiffusionDirection:BOOL=OFF
  -DANTS_BUILD_antsSliceRegularizedRegistration:BOOL=OFF
  -DANTS_BUILD_ANTSIntegrateVectorField:BOOL=OFF
  -DANTS_BUILD_ANTSIntegrateVelocityField:BOOL=OFF
  -DANTS_BUILD_antsTransformInfo:BOOL=OFF
  -DANTS_BUILD_antsUtilitiesTesting:BOOL=OFF
  -DANTS_BUILD_AverageAffineTransform:BOOL=OFF
  -DANTS_BUILD_AverageAffineTransformNoRigid:BOOL=OFF
  -DANTS_BUILD_AverageImages:BOOL=OFF
  -DANTS_BUILD_simpleSynRegistration:BOOL=OFF
  -DANTS_BUILD_CompositeTransformUtil:BOOL=OFF
  -DANTS_BUILD_CreateDisplacementField:BOOL=OFF
  -DANTS_BUILD_ConvertTransformFile:BOOL=OFF
  -DANTS_BUILD_compareTwoTransforms:BOOL=OFF
  -DANTS_BUILD_SuperResolution:BOOL=OFF
  -DANTS_BUILD_TimeSCCAN:BOOL=OFF
  -DANTS_BUILD_TextureCooccurrenceFeatures:BOOL=OFF
  -DANTS_BUILD_TextureRunLengthFeatures:BOOL=OFF
  -DANTS_BUILD_ImageIntensityStatistics:BOOL=OFF
  -DANTS_BUILD_GetConnectedComponentsFeatureImages:BOOL=OFF
  -DANTS_BUILD_DeNrrd:BOOL=OFF
  -DANTS_BUILD_StackSlices:BOOL=OFF
  )
if(${SUPERBUILD_TOPLEVEL_PROJECT}_USE_QT)
  list(APPEND ${proj}_CMAKE_OPTIONS -DANTS_USE_QT:BOOL=ON)
endif()
### --- End Project specific additions
set(${proj}_REPOSITORY "https://github.com/ANTsX/ANTs.git")
set(${proj}_GIT_TAG
  51a52c994b22081a68f4102301e034ec45cfd4f7 # 20210609
)

ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  GIT_REPOSITORY ${${proj}_REPOSITORY}
  GIT_TAG ${${proj}_GIT_TAG}
  SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
  BINARY_DIR ${proj}-${EXTERNAL_PROJECT_BUILD_TYPE}-build
  LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
  LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
  LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
  LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
  ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS -Wno-dev --no-warn-unused-cli
  CMAKE_CACHE_ARGS
    ${${proj}_CMAKE_OPTIONS}
    ${EXTERNAL_PROJECT_DEFAULTS}
  #INSTALL_COMMAND ""
  DEPENDS "${${proj}_DEPENDENCIES}"
  )

set(${proj}_SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj})
set(${proj}_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/lib)
#${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}-${CMAKE_BUILD_TYPE}-EP${EXTERNAL_PROJECT_BUILD_TYPE}-build/lib)

mark_as_superbuild(
  VARS
     ${proj}_SOURCE_DIR:PATH
     ${proj}_LIBRARY_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
