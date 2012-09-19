
add_executable(PrettyPrintTableTest PrettyPrintTableTest.cxx)
set_target_properties(PrettyPrintTableTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin)

add_test(NAME PrettyPrintTableTest COMMAND ${LAUNCH_EXE} $<TARGET_FILE:PrettyPrintTableTest>)

add_executable(BRAINSCleanMask BRAINSCleanMask.cxx)
target_link_libraries(BRAINSCleanMask ${ITK_LIBRARIES})

## HACK: TODO: Make this test do something meaningful.
#add_executable( itkResampleInPlaceImageFilterTest itkBasicFiltersTests.cxx )
#target_link_libraries( itkResampleInPlaceImageFilterTest ITKIO ITKBasicFilters )
#add_test( itkResampleInPlaceImageFilterTest itkResampleInPlaceImageFilterTest
#    itkResampleInPlaceImageFilterTest
#    ${ITK_DATA_ROOT}/Input/itkResampleInPlaceImageFilterTest.nii.gz
#    ${ITK_DATA_ROOT}/Baseline/Review/itkResampleInPlaceImageFilterTest.nii.gz
#    ${CMAKE_CURRENT_BINARY_DIR}/test_image.nii.gz
#)

