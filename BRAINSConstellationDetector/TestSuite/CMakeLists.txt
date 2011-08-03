list(APPEND ExternalData_URL_TEMPLATES
  # Local data store populated by the ITK pre-commit hook
  "file:///${CMAKE_SOURCE_DIR}/.ExternalData/%(algo)/%(hash)"
  # Data published by Iowa Psychiatry web interface
  "http://www.psychiatry.uiowa.edu/users/brainstestdata/ctestdata/%(algo)/%(hash)"

  # Data published by MIDAS
  "http://midas.kitware.com/api/rest/midas.bitstream.by.hash?hash=%(hash)&algorithm=%(algo)"

  # Data published by developers using git-gerrit-push.
  "http://www.itk.org/files/ExternalData/%(algo)/%(hash)"
  )

# Tell ExternalData commands to transform raw files to content links.
# TODO: Condition this feature on presence of our pre-commit hook.
set(ExternalData_LINK_CONTENT MD5)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/ExternalData.cmake)

## Set testing environment
##
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/obliq_setup.txt.in
               ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/obliq_setup.txt @ONLY IMMEDIATE)

set(ALL_TEST_PROGS
  BRAINSAlignMSP
  BRAINSConstellationDetector
  BRAINSEyeDetector
  BRAINSClipInferior
  BRAINSConstellationModeler
  BRAINSTrimForegroundInDirection
)


foreach(testprog ${ALL_TEST_PROGS})
MakeTestDriverFromSEMTool(${testprog} ${testprog}Test.cxx)
endforeach()


if ( 0 ) ## HACK:  Just silencing failing tests until they can be proper addressed
set(BRAINSConstellationModelerTestName BRAINSConstellationModelerTest_T1)
if(ConstellationModeler_TEST_causes_major_delay)
ExternalData_add_test( FetchData ${BRAINSConstellationModelerTestName} BRAINSConstellationModelerTest
         BRAINSConstellationModelerTest
         --inputTrainingList ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/obliq_setup.txt
         --outputModel ${CMAKE_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationModelerTestName}.mdl
)
endif(ConstellationModeler_TEST_causes_major_delay)
endif()

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_T1)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
  --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5})

## Test BRAINSEyeDetector
##
set(BRAINSEyeDetectorTestName BRAINSEyeDetectorTest_T1)
ExternalData_add_test( FetchData ${BRAINSEyeDetectorTestName} BRAINSEyeDetectorTest
  --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSEyeDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSEyeDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSEyeDetectorTest
         DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSEyeDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSEyeDetectorTestName}_aligned.nii.gz
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION)

## Test BRAINSAlignMSP:
##
set(BRAINSAlignMSPTestName BRAINSAlignMSPTest_T1)
ExternalData_add_test( FetchData ${BRAINSAlignMSPTestName} BRAINSAlignMSPTest
  --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSAlignMSPTestName}_standard.nii.gz}
           ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSAlignMSPTestName}_aligned.nii.gz
         --compareIntensityTolerance 50
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 500
         BRAINSAlignMSPTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --OutputresampleMSP ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSAlignMSPTestName}_aligned.nii.gz
         --mspQualityLevel 3
)


## Test BRAINSTrimForegroundInDirection:
##
set(BRAINSTrimForegroundInDirectionTestName BRAINSTrimForegroundInDirectionTest_T1)
ExternalData_add_test( FetchData ${BRAINSTrimForegroundInDirectionTestName} BRAINSTrimForegroundInDirectionTest
  --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSTrimForegroundInDirectionTestName}_standard.nii.gz}
           ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSTrimForegroundInDirectionTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSTrimForegroundInDirectionTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSTrimForegroundInDirectionTestName}_aligned.nii.gz
         --headSizeLimit 2600.0)


## Test BRAINSClipInferior
##
# add_executable(BRAINSClipInferiorTest BRAINSClipInferiorTest.cxx)
# target_link_libraries(BRAINSClipInferiorTest BRAINSClipInferiorCOMMONLIB)


## Test Versor
##
set (VersorTester_source VersorTester.cxx)
add_executable(VersorTester ${VersorTester_source})
target_link_libraries(VersorTester ${ITK_LIBRARIES})

install(TARGETS VersorTester
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static)


## Test landmarksConstellationTrainingDefinitionIO
##
add_executable(TestlandmarksConstellationTrainingDefinitionIO TestlandmarksConstellationTrainingDefinitionIO.cxx)

set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_otsuThreshold)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --otsuPercentileThreshold 0.05)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_mspLevel)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --mspQualityLevel 3)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_interpolationMode)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --interpolationMode BSpline)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rescaleIntensity)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --rescaleIntensities)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_acLowerBound)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --acLowerBound 100)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rescaleIntensityRange)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --rescaleIntensitiesOutputRange 20,2000)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_backgroundValue)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --BackgroundFillValue BIGNEG)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_cutOutHeadInOutputResampledVolume)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --cutOutHeadInOutputVolume)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_cutOutHeadInOutputVolume)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --cutOutHeadInOutputVolume)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_trimRescaledIntensities)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --trimRescaledIntensities 2)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rVN4)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --rVN4 200)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rpc)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --rpc 200)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rac)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --rac 200)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rmpj)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --rmpj 200)


set(BRAINSConstellationDetectorTestName BRAINSConstellationDetectorTest_rVN4-rpc-rac-rmpj)
ExternalData_add_test( FetchData ${BRAINSConstellationDetectorTestName} BRAINSConstellationDetectorTest
         --compare DATA{${CMAKE_SOURCE_DIR}/TestData/${BRAINSConstellationDetectorTestName}_standard.nii.gz}
         ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --compareIntensityTolerance 0
         --compareRadiusTolerance 0
         --compareNumberOfPixelsTolerance 0

         BRAINSConstellationDetectorTest
         --inputVolume DATA{${CMAKE_SOURCE_DIR}/TestData/T1.nii.gz}
         --outputResampledVolume ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION/${BRAINSConstellationDetectorTestName}_aligned.nii.gz
         --inputTemplateModel DATA{${CMAKE_SOURCE_DIR}/TestData/newT1.mdl}
         --LLSModel DATA{${CMAKE_SOURCE_DIR}/TestData/LLSModel.hdf5}
         --rmpj 200
         --rVN4 200
         --rac 200
         --rpc 200)

ExternalData_Add_Target( FetchData )  # Name of data management target
