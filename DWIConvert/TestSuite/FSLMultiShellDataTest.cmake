# FSL Multishell data test
# first, convert multishell FSL into Nrrd
set(command_line
        ${TEST_PROGRAM}
        --inputVolume ${INPUTVOLUME}
        --outputVolume ${OUTPUTVOLUME}
        )

message("Running ${command_line}")

execute_process(COMMAND ${command_line}
        RESULT_VARIABLE TEST_RESULT
        )

if(TEST_RESULT)
    message(FATAL_ERROR "Converting multishell FSL into Nrrd failed")
endif()

# second, convert Nrrd back to FSL
set(command_line
        ${TEST_PROGRAM}
        --inputVolume ${OUTPUTVOLUME}
        --outputVolume ${RECOVERVOLUME}
        )

message("Running ${command_line}")

execute_process(COMMAND ${command_line}
        RESULT_VARIABLE TEST_RESULT
        )

if(TEST_RESULT)
    message(FATAL_ERROR "Recovering  multishell Nrrd back to FSL failed")
endif()

# compare the bvector and bvalue of recovering FSL files with original FSL files
set(command_line
        ${TEXT_COMPARE_PROGRAM}
        ${INPUTBVECTOR}
        ${RECOVERBVECTOR}
        ${INPUTBVALUE}
        ${RECOVERBVALUE}
        )
message("Running ${command_line}")
execute_process(COMMAND
        ${command_line}
        RESULT_VARIABLE TEST_RESULT
        )
if(TEST_RESULT)
    message(FATAL_ERROR
            "Failed: Recovering bvec or bval files do not match the original FSL file")
endif()

message("Passed")