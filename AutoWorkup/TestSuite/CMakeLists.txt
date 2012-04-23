## Set testing environment
##

## TODO:  This needs to become a real test where several configure file options are used
##        to transform the hardcoded values in the scripts in this directory to build
##        specific versions.
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME AutoWorkupDownloadData
    COMMAND echo
    GetAutoWorkupDataTest
    DATA{${TestData_DIR}/SUBJ_A_small_T1.nii.gz}
    DATA{${TestData_DIR}/SUBJ_A_small_T2.nii.gz}
    DATA{${TestData_DIR}/SUBJ_B_small_T1.nii.gz}
    DATA{${TestData_DIR}/SUBJ_B_small_T2.nii.gz}
)

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
