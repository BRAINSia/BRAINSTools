#
# first reconstruct a nrrd file from the generated FSL files
set(command_line
  ${FSL_TO_NRRD}
  --conversionMode FSLToNrrd
  --inputBVectors ${VEC_FILE}
  --inputBValues ${VAL_FILE}
  --inputVolume ${NII_FILE}
  --outputVolume ${NRRD_FILE}
)

message("Running ${command_line}")

execute_process(COMMAND ${command_line}
  RESULT_VARIABLE TEST_RESULT
)

if(TEST_RESULT)
  message(FATAL_ERROR "${FSL_TO_NRRD} run failed")
endif()

#
# compare the results

set(command_line
  ${TEST_COMPARE_PROGRAM} --inputVolume1 ${NRRD_FILE}
  --inputVolume2 ${NRRD_COMPARE_FILE}
)

message("Running ${command_line}")

execute_process(COMMAND
  ${command_line}
  RESULT_VARIABLE TEST_RESULT
)

if(TEST_RESULT)
  message(FATAL_ERROR
    "Failed: ${NRRD_FILE} doesn't match ${NRRD_COMPARE_FILE}")
endif()

message("Passed")
