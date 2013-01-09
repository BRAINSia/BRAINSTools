set(IMAGE_PATH_OUTPUTS ${ConvertBetweenFileFormats_BINARY_DIR})
# Convert from grayscale mhd to other image file formats
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_to_png COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_out.png)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_to_jpg_unsigned_int COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_out.jpg unsigned_int)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_to_tiff COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_out.tiff)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_to_analyze75 COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> DATA{${TestData_DIR}/BigEndian.mhd,BigEndian.raw} ${IMAGE_PATH_OUTPUTS}/image_out.img)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_to_bmp_unsigned_char COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_out.bmp unsigned_char)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_to_nrrd_float COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_out.nrrd float)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_to_vtk_long COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_out.vtk long)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_to_gipl COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_out.gipl)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_to_nifti COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_out.nii.gz)

# Convert back to mhd
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_from_png COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> ${IMAGE_PATH_OUTPUTS}/image_out.png    ${IMAGE_PATH_OUTPUTS}/image_in1.mhd)
set_tests_properties(ccvnt_from_png PROPERTIES DEPENDS ccvnt_to_png)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_from_jpg_unsigned_int COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> ${IMAGE_PATH_OUTPUTS}/image_out.jpg    ${IMAGE_PATH_OUTPUTS}/image_in2.mhd unsigned_short)
set_tests_properties(ccvnt_from_jpg_unsigned_int PROPERTIES DEPENDS ccvnt_to_jpg_unsigned_int)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_from_tiff COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> ${IMAGE_PATH_OUTPUTS}/image_out.tiff   ${IMAGE_PATH_OUTPUTS}/image_in3.mhd)
set_tests_properties(ccvnt_from_tiff PROPERTIES DEPENDS ccvnt_to_tiff)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_from_analyze75 COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> ${IMAGE_PATH_OUTPUTS}/image_out.hdr    ${IMAGE_PATH_OUTPUTS}/image_in4.mhd)
set_tests_properties(ccvnt_from_analyze75 PROPERTIES DEPENDS ccvnt_to_analyze75)

#
# Disabling this test because BMP are now only supporting reading of RGB images.
# ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_from_bmp_unsigned_char COMMAND
#          $<TARGET_FILE:ConvertBetweenFileFormats> ${IMAGE_PATH_OUTPUTS}/image_out.bmp    ${IMAGE_PATH_OUTPUTS}/image_in5.mhd unsigned_short)
#
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_from_nrrd_float COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> ${IMAGE_PATH_OUTPUTS}/image_out.nrrd   ${IMAGE_PATH_OUTPUTS}/image_in6.mhd unsigned_short)
set_tests_properties(ccvnt_from_nrrd_float PROPERTIES DEPENDS ccvnt_to_nrrd_float)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_from_vtk_long COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> ${IMAGE_PATH_OUTPUTS}/image_out.vtk    ${IMAGE_PATH_OUTPUTS}/image_in7.mhd unsigned_short)
set_tests_properties(ccvnt_from_vtk_long PROPERTIES DEPENDS ccvnt_to_vtk_long)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_from_gipl COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> ${IMAGE_PATH_OUTPUTS}/image_out.gipl   ${IMAGE_PATH_OUTPUTS}/image_in8.mhd)
set_tests_properties(ccvnt_from_gipl PROPERTIES DEPENDS ccvnt_to_gipl)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ccvnt_from_nifti COMMAND
          $<TARGET_FILE:ConvertBetweenFileFormats> ${IMAGE_PATH_OUTPUTS}/image_out.nii.gz ${IMAGE_PATH_OUTPUTS}/image_in8.mhd)
set_tests_properties(ccvnt_from_nifti PROPERTIES DEPENDS ccvnt_to_nifti)

# Compare with original png
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME compare_png COMMAND
          $<TARGET_FILE:ImageCompareTests> --compare DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_in1.mhd Dummy)
set_tests_properties(compare_png PROPERTIES DEPENDS ccvnt_from_png)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME compare_jpg_unsigned_int COMMAND
          $<TARGET_FILE:ImageCompareTests> DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_in2.mhd)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME compare_tiff COMMAND
          $<TARGET_FILE:ImageCompareTests> --compare DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_in3.mhd Dummy)
set_tests_properties(compare_tiff PROPERTIES DEPENDS ccvnt_from_tiff)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME compare_analyze75 COMMAND
          $<TARGET_FILE:ImageCompareTests> --compare DATA{${TestData_DIR}/BigEndian.mhd,BigEndian.raw} ${IMAGE_PATH_OUTPUTS}/image_in4.mhd Dummy)
set_tests_properties(compare_analyze75 PROPERTIES DEPENDS ccvnt_from_analyze75)

#
# Disabling this test because BMP are now only supporting reading of RGB images.
# ExternalData_add_test( ${PROJECT_NAME}FetchData NAME compare_bmp_unsigned_char COMMAND
#          $<TARGET_FILE:ImageCompareTests> --compare DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_in5.mhd Dummy)
#
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME compare_nrrd_float COMMAND
          $<TARGET_FILE:ImageCompareTests> --compare DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_in6.mhd Dummy)
set_tests_properties(compare_nrrd_float PROPERTIES DEPENDS ccvnt_from_nrrd_float)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME compare_vtk_long COMMAND
          $<TARGET_FILE:ImageCompareTests> --compare DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_in7.mhd Dummy)
set_tests_properties(compare_vtk_long PROPERTIES DEPENDS ccvnt_from_vtk_long)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME compare_gipl COMMAND
          $<TARGET_FILE:ImageCompareTests> --compare DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_in8.mhd Dummy)
set_tests_properties(compare_gipl PROPERTIES DEPENDS ccvnt_from_gipl)

ExternalData_add_test( ${PROJECT_NAME}FetchData NAME compare_nifti COMMAND
          $<TARGET_FILE:ImageCompareTests> --compare DATA{${TestData_DIR}/image_in.mhd,image_in.raw} ${IMAGE_PATH_OUTPUTS}/image_in8.mhd Dummy)
set_tests_properties(compare_nifti PROPERTIES DEPENDS ccvnt_from_nifti)

add_executable(ImageCompareTests ImageCompareTests.cxx)
target_link_libraries(ImageCompareTests ${ITK_LIBRARIES})

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target

