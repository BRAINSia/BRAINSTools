MakeTestDriverFromSEMTool(BRAINSMultiSTAPLE BRAINSMultiSTAPLETest.cxx)

#if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 5) #This test takes way to long to run at the current scale.  It needs to be made smaller
#ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET}
#  NAME BRAINSMultiSTAPLETest
#  COMMAND $<TARGET_FILE:BRAINSMultiSTAPLE>
#  --inputLabelVolume DATA{${TestData_DIR}/BRAINSMultiSTAPLE/LM_1_seg.nii.gz}
#  --inputLabelVolume DATA{${TestData_DIR}/BRAINSMultiSTAPLE/LM_2_seg.nii.gz}
#  --inputLabelVolume DATA{${TestData_DIR}/BRAINSMultiSTAPLE/LM_3_seg.nii.gz}
#  --inputTransform DATA{${TestData_DIR}/BRAINSMultiSTAPLE/session1composite_1.h5}
#  --inputTransform DATA{${TestData_DIR}/BRAINSMultiSTAPLE/session1composite_2.h5}
#  --inputTransform DATA{${TestData_DIR}/BRAINSMultiSTAPLE/session1composite_3.h5}
#  --inputCompositeT1Volume DATA{${TestData_DIR}/BRAINSMultiSTAPLE/AVG_T1.nii.gz}
#  --outputConfusionMatrix ${CMAKE_CURRENT_BINARY_DIR}/Confusion.mat
#  --outputMultiSTAPLE ${CMAKE_CURRENT_BINARY_DIR}/MultiSTAPLE.nii.gz
#  )
#endif()
