MakeTestDriverFromSEMTool(BRAINSMultiSTAPLE
  BRAINSMultiSTAPLETest.cxx)

ExternalData_add_test(BRAINSMultiSTAPLEFetchData
  NAME BRAINSMultiSTAPLETest
  COMMAND $<TARGET_FILE:BRAINSMultiSTAPLE>
  --inputLabelVolume DATA{${TestData_DIR}/BRAINSMultiSTAPLE/LM_1_seg.nii.gz}
  --inputLabelVolume DATA{${TestData_DIR}/BRAINSMultiSTAPLE/LM_2_seg.nii.gz}
  --inputLabelVolume DATA{${TestData_DIR}/BRAINSMultiSTAPLE/LM_3_seg.nii.gz}
  --inputTransform DATA{${TestData_DIR}/BRAINSMultiSTAPLE/session1composite_1.h5}
  --inputTransform DATA{${TestData_DIR}/BRAINSMultiSTAPLE/session1composite_2.h5}
  --inputTransform DATA{${TestData_DIR}/BRAINSMultiSTAPLE/session1composite_3.h5}
  --inputCompositeT1Volume DATA{${TestData_DIR}/BRAINSMultiSTAPLE/AVG_T1.nii.gz}
  --outputConfusionMatrix ${CMAKE_CURRENT_BINARY_DIR}/Confusion.mat
  --outputMultiSTAPLE ${CMAKE_CURRENT_BINARY_DIR}/MultiSTAPLE.nii.gz
  )
ExternalData_Add_Target(BRAINSMultiSTAPLEFetchData)
