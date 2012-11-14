# This project is intended to be built outside the Insight source tree
PROJECT(ConvertBetweenFileFormats)

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/castconvertConfigure.h.in 
                ${CMAKE_CURRENT_BINARY_DIR}/castconvertConfigure.h @ONLY IMMEDIATE )

ADD_EXECUTABLE(ConvertBetweenFileFormats castconvert.cxx  castconvertDicomScalar.cxx castconvertScalar.cxx  castconvertScalar2D.cxx castconvertScalar3D.cxx castconvertScalar4D.cxx castconvertScalar4DA.cxx castconvertScalar2DA.cxx castconvertScalar3DA.cxx castconvertDicomScalarA.cxx)

INCLUDE_DIRECTORIES( 
    ${CMAKE_CURRENT_SOURCE_DIR} 
    ${CMAKE_CURRENT_BINARY_DIR} 
    )

TARGET_LINK_LIBRARIES(ConvertBetweenFileFormats  ${ITK_LIBRARIES} )
IF (VTK_FOUND)
  TARGET_LINK_LIBRARIES(ConvertBetweenFileFormats vtkIO vtkCommon vtkFiltering vtkImaging )
  INCLUDE_DIRECTORIES( ${ITKApps_SOURCE_DIR}/Auxiliary/vtk )
ENDIF (VTK_FOUND)

INSTALL_TARGETS(/bin ConvertBetweenFileFormats)

IF( BUILD_TESTING )
  SET(CXX_TEST_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  SET(CONVERT_TESTS ${CXX_TEST_PATH}/ConvertBetweenFileFormats)

  SET(IMAGE_PATH ${ConvertBetweenFileFormats_SOURCE_DIR}/Testing/Data)
  SET(IMAGE_PATH_OUTPUTS ${ConvertBetweenFileFormats_BINARY_DIR})
  # Convert from grayscale mhd to other image file formats
  ADD_TEST(ccvnt_to_png                   ${CONVERT_TESTS} ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_out.png)
  ADD_TEST(ccvnt_to_jpg_unsigned_int      ${CONVERT_TESTS} ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_out.jpg unsigned_int)
  ADD_TEST(ccvnt_to_tiff                  ${CONVERT_TESTS} ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_out.tiff)
  ADD_TEST(ccvnt_to_analyze75             ${CONVERT_TESTS} ${IMAGE_PATH}/BigEndian.mhd ${IMAGE_PATH_OUTPUTS}/image_out.img)
  ADD_TEST(ccvnt_to_bmp_unsigned_char     ${CONVERT_TESTS} ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_out.bmp unsigned_char)
  ADD_TEST(ccvnt_to_nrrd_float            ${CONVERT_TESTS} ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_out.nrrd float)
  ADD_TEST(ccvnt_to_vtk_long              ${CONVERT_TESTS} ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_out.vtk long)
  ADD_TEST(ccvnt_to_gipl                  ${CONVERT_TESTS} ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_out.gipl)
  ADD_TEST(ccvnt_to_nifti                 ${CONVERT_TESTS} ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_out.nii.gz)

  # Convert back to mhd
  ADD_TEST(ccvnt_from_png                 ${CONVERT_TESTS} ${IMAGE_PATH_OUTPUTS}/image_out.png    ${IMAGE_PATH_OUTPUTS}/image_in1.mhd)
  set_tests_properties(ccvnt_from_png PROPERTIES DEPENDS ccvnt_to_png)

  ADD_TEST(ccvnt_from_jpg_unsigned_int    ${CONVERT_TESTS} ${IMAGE_PATH_OUTPUTS}/image_out.jpg    ${IMAGE_PATH_OUTPUTS}/image_in2.mhd unsigned_short)
  set_tests_properties(ccvnt_from_jpg_unsigned_int PROPERTIES DEPENDS ccvnt_to_jpg_unsigned_int)

  ADD_TEST(ccvnt_from_tiff                ${CONVERT_TESTS} ${IMAGE_PATH_OUTPUTS}/image_out.tiff   ${IMAGE_PATH_OUTPUTS}/image_in3.mhd)
  set_tests_properties(ccvnt_from_tiff PROPERTIES DEPENDS ccvnt_to_tiff)

  ADD_TEST(ccvnt_from_analyze75           ${CONVERT_TESTS} ${IMAGE_PATH_OUTPUTS}/image_out.hdr    ${IMAGE_PATH_OUTPUTS}/image_in4.mhd)
  set_tests_properties(ccvnt_from_analyze75 PROPERTIES DEPENDS ccvnt_to_analyze75)

  #
  # Disabling this test because BMP are now only supporting reading of RGB images.
  # ADD_TEST(ccvnt_from_bmp_unsigned_char   ${CONVERT_TESTS} ${IMAGE_PATH_OUTPUTS}/image_out.bmp    ${IMAGE_PATH_OUTPUTS}/image_in5.mhd unsigned_short)
  #
  ADD_TEST(ccvnt_from_nrrd_float          ${CONVERT_TESTS} ${IMAGE_PATH_OUTPUTS}/image_out.nrrd   ${IMAGE_PATH_OUTPUTS}/image_in6.mhd unsigned_short)
  set_tests_properties(ccvnt_from_nrrd_float PROPERTIES DEPENDS ccvnt_to_nrrd_float)

  ADD_TEST(ccvnt_from_vtk_long            ${CONVERT_TESTS} ${IMAGE_PATH_OUTPUTS}/image_out.vtk    ${IMAGE_PATH_OUTPUTS}/image_in7.mhd unsigned_short)
  set_tests_properties(ccvnt_from_vtk_long PROPERTIES DEPENDS ccvnt_to_vtk_long)

  ADD_TEST(ccvnt_from_gipl                ${CONVERT_TESTS} ${IMAGE_PATH_OUTPUTS}/image_out.gipl   ${IMAGE_PATH_OUTPUTS}/image_in8.mhd)
  set_tests_properties(ccvnt_from_gipl PROPERTIES DEPENDS ccvnt_to_gipl)

  ADD_TEST(ccvnt_from_nifti               ${CONVERT_TESTS} ${IMAGE_PATH_OUTPUTS}/image_out.nii.gz ${IMAGE_PATH_OUTPUTS}/image_in8.mhd)
  set_tests_properties(ccvnt_from_nifti PROPERTIES DEPENDS ccvnt_to_nifti)

  SET(IMAGE_COMPARE_TESTS ${CXX_TEST_PATH}/ImageCompareTests)

  # Compare with original png
  ADD_TEST(compare_png                    ${IMAGE_COMPARE_TESTS} --compare ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in1.mhd Dummy)
  set_tests_properties(compare_png PROPERTIES DEPENDS ccvnt_from_png)

 #ADD_TEST(compare_jpg_unsigned_int       ${IMAGE_COMPARE_TESTS} ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in2.mhd)
  ADD_TEST(compare_tiff                   ${IMAGE_COMPARE_TESTS} --compare ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in3.mhd Dummy)
  set_tests_properties(compare_tiff PROPERTIES DEPENDS ccvnt_from_tiff)

  ADD_TEST(compare_analyze75              ${IMAGE_COMPARE_TESTS} --compare ${IMAGE_PATH}/BigEndian.mhd ${IMAGE_PATH_OUTPUTS}/image_in4.mhd Dummy)
  set_tests_properties(compare_analyze75 PROPERTIES DEPENDS ccvnt_from_analyze75)

  #
  # Disabling this test because BMP are now only supporting reading of RGB images.
  # ADD_TEST(compare_bmp_unsigned_char      ${IMAGE_COMPARE_TESTS} --compare ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in5.mhd Dummy)
  #
  ADD_TEST(compare_nrrd_float             ${IMAGE_COMPARE_TESTS} --compare ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in6.mhd Dummy)
  set_tests_properties(compare_nrrd_float PROPERTIES DEPENDS ccvnt_from_nrrd_float)

  ADD_TEST(compare_vtk_long               ${IMAGE_COMPARE_TESTS} --compare ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in7.mhd Dummy)
  set_tests_properties(compare_vtk_long PROPERTIES DEPENDS ccvnt_from_vtk_long)

  ADD_TEST(compare_gipl                   ${IMAGE_COMPARE_TESTS} --compare ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in8.mhd Dummy)
  set_tests_properties(compare_gipl PROPERTIES DEPENDS ccvnt_from_gipl)

  ADD_TEST(compare_nifti            ${IMAGE_COMPARE_TESTS} --compare ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in8.mhd Dummy)
  set_tests_properties(compare_nifti PROPERTIES DEPENDS ccvnt_from_nifti)

  ADD_EXECUTABLE(ImageCompareTests ImageCompareTests.cxx)
  TARGET_LINK_LIBRARIES(ImageCompareTests ${ITK_LIBRARIES})
ENDIF( BUILD_TESTING )

