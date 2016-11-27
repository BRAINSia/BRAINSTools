# arguments checking
if( NOT TEST_PROGRAM )
  message( FATAL_ERROR "Require TEST_PROGRAM to be defined" )
endif( )

if( NOT TEST_COMPARE_PROGRAM )
  message( FATAL_ERROR "Require TEST_COMPARE_PROGRAM to be defined" )
endif( )

if( NOT TEST_BASELINE )
  message( FATAL_ERROR "Require TEST_BASELINE to be defined" )
endif( )

if( NOT TEST_BASELINE_NHDR )
  set( TEST_BASELINE_NHDR ${TEST_BASELINE} )
endif( )

if( NOT TEST_INPUT )
  message( FATAL_ERROR "Require TEST_INPUT to be defined" )
endif( )

if( NOT TEST_TEMP_OUTPUT )
  message( FATAL_ERROR "Require TEST_TEMP_OUTPUT to be defined" )
endif( )

# Run the compare program to make sure it built correctly
# execute_process(
#   COMMAND ${TEST_COMPARE_PROGRAM} --help
#   RESULT_VARIABLE TEST_RESULT
#   )

# if the return value is !=0 bail out
if( TEST_RESULT )
  message( FATAL_ERROR "Failed: Test compare program ${TEST_COMPARE_PROGRAM} won't run.\n${TEST_ERROR}" )
endif( )

# Check to see if the image we are comparing against exists.  We do this here to avoid a lengthy test for no reason.
if(NOT EXISTS ${TEST_BASELINE})
  message( FATAL_ERROR "Failed: Baseline image ${TEST_BASELINE} does not exist!\n")
endif( )


# run the test program, capture the stdout/stderr and the result var
if(TEST_FSL_FLAG)
  set(Test_Command_Line
    ${TEST_PROGRAM} --inputDicomDirectory ${TEST_INPUT} --outputVolume ${TEST_TEMP_OUTPUT} ${TEST_PROGRAM_ARGS}
  )
  list(INSERT Test_Command_Line 1 --conversionMode)
  list(INSERT Test_Command_Line 2 DicomToFSL)

  message("Test_Command_Line=${Test_Command_Line}")
  execute_process(
    COMMAND ${Test_Command_Line}
    ERROR_VARIABLE TEST_ERROR
    RESULT_VARIABLE TEST_RESULT
  )
  if( TEST_RESULT )
    message( FATAL_ERROR "Failed: Test program ${TEST_PROGRAM} exited != 0.\n${TEST_ERROR}" )
  endif( )

else()
  ## For testing purposes, convert from dicom to nhdr, nhdr to fsl, fsl to nrrd, then compare results
  get_filename_component(TEST_TEMP_DIRNAME ${TEST_TEMP_OUTPUT} DIRECTORY)
  get_filename_component(TEST_TEMP_BASENAME ${TEST_TEMP_OUTPUT} NAME_WE )
  set(TEST_TEMP_FULLBASENAME "${TEST_TEMP_DIRNAME}/${TEST_TEMP_BASENAME}")


  set(TEST_INITIAL_FSL_OUTPUT ${TEST_TEMP_FULLBASENAME}_NRRDTOFSL.nii.gz)
  set(TEST_INITIAL_NHDR_OUTPUT ${TEST_TEMP_FULLBASENAME}_FSLTONHDR.nhdr)
  file(REMOVE ${TEST_INITIAL_NHDR_OUTPUT})
  file(REMOVE ${TEST_INITIAL_FSL_OUTPUT})
  file(REMOVE ${TEST_TEMP_OUTPUT})


  set(Test_Command_Line_DICOMToNRRD
    ${TEST_PROGRAM} --conversionMode DicomToNrrd --inputDicomDirectory ${TEST_INPUT} --outputVolume ${TEST_TEMP_OUTPUT} ${TEST_PROGRAM_ARGS}
  )
  set(Test_Command_Line_NRRDToFSL
    ${TEST_PROGRAM} --conversionMode NrrdToFSL --inputVolume ${TEST_TEMP_OUTPUT} --outputVolume ${TEST_INITIAL_FSL_OUTPUT} ${TEST_PROGRAM_ARGS}
  )
  set(Test_Command_Line_FSLToNHDR
    ${TEST_PROGRAM} --conversionMode FSLToNrrd --inputVolume ${TEST_INITIAL_FSL_OUTPUT} --outputVolume ${TEST_INITIAL_NHDR_OUTPUT} ${TEST_PROGRAM_ARGS}
  )

  message(STATUS "\n\nTest_Command_Line_DICOMToNRRD=${Test_Command_Line_DICOMToNRRD}\n")
  execute_process(
    COMMAND ${Test_Command_Line_DICOMToNRRD}
    ERROR_VARIABLE TEST_ERROR
    RESULT_VARIABLE TEST_RESULT
  )
  message(STATUS ${TEST_ERROR})
  if( TEST_RESULT )
    message( FATAL_ERROR "Failed: Test program ${Test_Command_Line_DICOMToNRRD} exited != 0.\n${TEST_ERROR}" )
  endif( )

  message(STATUS "\n\nTest_Command_Line_NRRDToFSL=${Test_Command_Line_NRRDToFSL}\n")
  execute_process(
    COMMAND ${Test_Command_Line_NRRDToFSL}
    ERROR_VARIABLE TEST_ERROR
    RESULT_VARIABLE TEST_RESULT
  )
  message(STATUS ${TEST_ERROR})
  if( TEST_RESULT )
    message( FATAL_ERROR "Failed: Test program ${Test_Command_Line_NRRDToFSL} exited != 0.\n${TEST_ERROR}" )
  endif( )

  message(STATUS "\n\nTest_Command_Line_FSLToNHDR=${Test_Command_Line_FSLToNHDR}\n")
  execute_process(
    COMMAND ${Test_Command_Line_FSLToNHDR}
    ERROR_VARIABLE TEST_ERROR
    RESULT_VARIABLE TEST_RESULT
  )
  message(STATUS ${TEST_ERROR})
  if( TEST_RESULT )
    message( FATAL_ERROR "Failed: Test program ${Test_Command_Line_FSLToNHDR} exited != 0.\n${TEST_ERROR}" )
  endif( )

endif()

if(TEST_FSL_FLAG)
    if(NOT TEST_BASELINE_BVEC)
      message( FATAL_ERROR "Require TEST_BASELINE_BVEC to be defined" )
    endif()
    if(NOT TEST_BASELINE_BVAL)
      message( FATAL_ERROR "Require TEST_BASELINE_BVAL to be defined" )
    endif()
    if(NOT TEST_TEMP_BVEC)
      message( FATAL_ERROR "Require TEST_TEMP_BVEC to be defined" )
    endif()
    if(NOT TEST_TEMP_BVAL)
      message( FATAL_ERROR "Require TEST_TEMP_BVAL to be defined" )
    endif()

    if(NOT TEST_TEXT_COMPARE)
      message( FATAL_ERROR "Require TEST_TEXT_COMPARE to be defined" )
    endif()
    set(Test_Compare_Command_Line_TEMP_OUTPUT
      ${TEST_TEXT_COMPARE} ${TEST_BASELINE_BVEC} ${TEST_TEMP_BVEC}
      ${TEST_BASELINE_BVAL} ${TEST_TEMP_BVAL})
    message("Test_Compare_Command_Line_TEMP_OUTPUT=${Test_Compare_Command_Line_TEMP_OUTPUT}")
    execute_process(COMMAND ${Test_Compare_Command_Line_TEMP_OUTPUT}
      ERROR_VARIABLE TEST_ERROR
      RESULT_VARIABLE TEST_RESULT)
    if( TEST_RESULT )
      message( FATAL_ERROR
        "Failed: Test program ${TEST_TEXT_COMPARE} exited != 0.  ${TEST_ERROR}" )
    endif( )
  #   foreach(file VEC VAL)
  #     file(READ ${TEST_BASELINE_B${file}} baseline)
  #     file(READ ${TEST_TEMP_B${file}} testoutput)
  #     if(testoutput STREQUAL baseline)
  #       message("B${file} output matches")
  #     else()
  #       message(FATAL_ERROR "B${file} output doesn't match
  # ${testoutput}
  # ${baseline}")
  #     endif()
  #   endforeach()
else()

  # if the return value is !=0 bail out
  if( TEST_RESULT )
    message( FATAL_ERROR "Failed: Test program ${TEST_PROGRAM} exited != 0.\n${TEST_ERROR}" )
  endif( )

  #------------------------------------------------------------------------------------------
  set(Test_Compare_Command_Line_TEMP_OUTPUT
    ${TEST_COMPARE_PROGRAM} --inputVolume2 ${TEST_TEMP_OUTPUT} --inputVolume1 ${TEST_BASELINE}
    )
  message("Test_Compare_Command_Line_TEMP_OUTPUT=${Test_Compare_Command_Line_TEMP_OUTPUT}")
  # now compare the output with the reference
  execute_process(
    COMMAND ${Test_Compare_Command_Line_TEMP_OUTPUT}
    RESULT_VARIABLE TEST_RESULT
    )
  # again, if return value is !=0 scream and shout
  if( TEST_RESULT )
    message( FATAL_ERROR "Failed: ${Test_Compare_Command_Line_TEMP_OUTPUT}: ${TEST_RESULT}")
  endif( )

  if( 0 ) ## These NHDR images are oriented differently, and would need a new BASELINE_NHDR ##TODO:
  #------------------------------------------------------------------------------------------
  set(Test_Compare_Command_Line_NHDR_OUTPUT
    ${TEST_COMPARE_PROGRAM} --inputVolume2 ${TEST_INITIAL_NHDR_OUTPUT} --inputVolume1 ${TEST_BASELINE_NHDR}
    )
  message("Test_Compare_Command_Line_NHDR_OUTPUT=${Test_Compare_Command_Line_NHDR_OUTPUT}")
  # now compare the output with the reference
  execute_process(
    COMMAND ${Test_Compare_Command_Line_NHDR_OUTPUT}
    RESULT_VARIABLE TEST_RESULT
    )
  # again, if return value is !=0 scream and shout
  if( TEST_RESULT )
    message( FATAL_ERROR "Failed: ${Test_Compare_Command_Line_NHDR_OUTPUT}: ${TEST_RESULT}")
  endif( )
  endif()
endif()

# everything went fine...
message( "Passed: The output of ${TEST_PROGRAM} matches ${TEST_BASELINE}:${TEST_BASELINE_NHDR}" )
