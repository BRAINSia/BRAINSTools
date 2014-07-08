#
# first, construct FSL files from NRRD
if(VEC_FILE AND VAL_FILE)
  set(command_line
    ${NRRD_TO_FSL}
    --conversionMode NrrdToFSL
    --inputVolume ${NRRD_FILE}
    --outputVolume ${NII_FILE}
    --outputBVectors ${VEC_FILE}
    --outputBValues ${VAL_FILE}
    )
else()
  # testing the case where the BVec/BVal filenames should
  # follow the output volume name
  set(command_line
    ${NRRD_TO_FSL}
    --conversionMode NrrdToFSL
    --inputVolume ${NRRD_FILE}
    --outputVolume ${NII_FILE}
    )
  string(REGEX REPLACE ".nii.gz$" ".bval" VAL_FILE "${NII_FILE}")
  string(REGEX REPLACE ".nii.gz$" ".bvec" VEC_FILE "${NII_FILE}")
endif()

message("Running ${command_line}")

execute_process(COMMAND ${command_line}
  RESULT_VARIABLE TEST_RESULT
)

if(TEST_RESULT)
  message(FATAL_ERROR "${NRRD_TO_FSL} run failed")
endif()

#
# compare NII files
set(command_line
  ${TEST_COMPARE_PROGRAM} --inputVolume1 ${NII_FILE}
  --inputVolume2 ${NII_COMPARE_FILE})

message("Running ${command_line}")

execute_process(COMMAND
  ${command_line}
  RESULT_VARIABLE TEST_RESULT
)

if(TEST_RESULT)
  message(FATAL_ERROR
    "Failed: ${NII_FILE} doesn't match ${NII_COMPARE_FILE}")
endif()

set(command_line
  ${TEXT_COMPARE_PROGRAM}
  ${VEC_FILE} ${VEC_COMPARE_FILE}
  ${VAL_FILE} ${VAL_COMPARE_FILE}
)
message("Running ${command_line}")

execute_process(COMMAND
  ${command_line}
  RESULT_VARIABLE TEST_RESULT
)

if(TEST_RESULT)
  message(FATAL_ERROR
    "Failed: FSL text files don't match")
endif()



message("Passed")
