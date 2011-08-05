project(BRAINSCUT_REGISTRATIONSCRIPTS)
##### Shell Scripts (TODO:: Eventually to the Python
configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/GenerateBSplineTransform_ROIAUTO.sh.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/GenerateBSplineTransform_ROIAUTO.sh IMMEDIATE @ONLY)
configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/GenerateBSplineTransform_ROI.sh.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/GenerateBSplineTransform_ROI.sh IMMEDIATE @ONLY)


configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/GenerateThDemon_CL_Transform_Gen_clipped_SumGrad_BSplineIniitilization.tcl.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/GenerateThDemon_CL_Transform_Gen_clipped_SumGrad_BSplineIniitilization.tcl IMMEDIATE @ONLY)

configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/GenerateBSplineTransform_ROI.tcl.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/GenerateBSplineTransform_ROI.tcl IMMEDIATE @ONLY)
configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/GenerateBSplineTransform_ROIAUTO.tcl.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/GenerateBSplineTransform_ROIAUTO.tcl IMMEDIATE @ONLY)
#configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/GenerateBSplineTransform_2010Jul.tcl.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/GenerateBSplineTransform_2010Jul.tcl IMMEDIATE @ONLY)

configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/GenerateDemonWarpTransform.tcl.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/GenerateDemonWarpTransform.tcl IMMEDIATE @ONLY)

configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/CheckOutputsNewer.tcl.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/CheckOutputsNewer.tcl IMMEDIATE @ONLY)

configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/GenerateThDemon_T1.tcl.in  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/GenerateThDemon_T1.tcl IMMEDIATE @ONLY)

configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/GenerateDemonWarpTransform_T1.tcl.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/GenerateDemonWarpTransform_T1.tcl IMMEDIATE @ONLY)

configure_file(${BRAINSCUT_REGISTRATIONSCRIPTS_SOURCE_DIR}/GenerateBSpline_BrainClipped.tcl.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/GenerateBSpline_BrainClipped.tcl IMMEDIATE @ONLY)
