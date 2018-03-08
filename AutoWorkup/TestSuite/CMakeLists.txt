## Set testing environment
##

if( ${BRAINSTools_MAX_TEST_LEVEL} GREATER 6)
## TODO:  This needs to become a real test where several configure file options are used
##        to transform the hardcoded values in the scripts in this directory to build
##        specific versions.
ExternalData_add_test( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} NAME AutoWorkupDownloadData
    COMMAND echo
    GetAutoWorkupDataTest
    DATA{${TestData_DIR}/SUBJ_A_small_T1.nii.gz}
    DATA{${TestData_DIR}/SUBJ_A_small_T2.nii.gz}
    DATA{${TestData_DIR}/SUBJ_B_small_T1.nii.gz}
    DATA{${TestData_DIR}/SUBJ_B_small_T2.nii.gz}
)
endif()
