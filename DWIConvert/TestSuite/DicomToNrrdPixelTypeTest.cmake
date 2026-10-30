#
# Convert a DICOM series to NRRD, round-trip it through FSL, and check that the
# NRRD "type:" field is EXPECTED_TYPE both times.
#
foreach(required TEST_PROGRAM TEST_INPUT TEST_TEMP_OUTPUT EXPECTED_TYPE)
  if(NOT ${required})
    message(FATAL_ERROR "Require ${required} to be defined")
  endif()
endforeach()

get_filename_component(TEST_TEMP_DIRNAME ${TEST_TEMP_OUTPUT} DIRECTORY)
get_filename_component(TEST_TEMP_BASENAME ${TEST_TEMP_OUTPUT} NAME_WE)
set(TEST_FSL_OUTPUT ${TEST_TEMP_DIRNAME}/${TEST_TEMP_BASENAME}_RoundTrip.nii.gz)
set(TEST_NHDR_OUTPUT ${TEST_TEMP_DIRNAME}/${TEST_TEMP_BASENAME}_RoundTrip.nhdr)
file(REMOVE ${TEST_TEMP_OUTPUT} ${TEST_FSL_OUTPUT} ${TEST_NHDR_OUTPUT})

function(run_dwiconvert)
  execute_process(
    COMMAND ${TEST_PROGRAM} ${ARGN} ${TEST_PROGRAM_ARGS}
    ERROR_VARIABLE TEST_ERROR
    RESULT_VARIABLE TEST_RESULT
  )
  if(TEST_RESULT)
    message(FATAL_ERROR "Failed: ${TEST_PROGRAM} ${ARGN} exited != 0.\n${TEST_ERROR}")
  endif()
endfunction()

function(check_nrrd_type nrrd_file)
  file(STRINGS ${nrrd_file} type_lines REGEX "^type:" LIMIT_COUNT 1)
  if(NOT type_lines STREQUAL "type: ${EXPECTED_TYPE}")
    message(FATAL_ERROR "${nrrd_file}: expected 'type: ${EXPECTED_TYPE}', found '${type_lines}'")
  endif()
endfunction()

run_dwiconvert(--conversionMode DicomToNrrd --inputDicomDirectory ${TEST_INPUT} --outputVolume ${TEST_TEMP_OUTPUT})
check_nrrd_type(${TEST_TEMP_OUTPUT})

run_dwiconvert(--conversionMode NrrdToFSL --inputVolume ${TEST_TEMP_OUTPUT} --outputVolume ${TEST_FSL_OUTPUT})
run_dwiconvert(--conversionMode FSLToNrrd --inputVolume ${TEST_FSL_OUTPUT} --outputVolume ${TEST_NHDR_OUTPUT})
check_nrrd_type(${TEST_NHDR_OUTPUT})
