# Dicom file has different zSpace and Thickness data test
# first, convert dicom data into Nrrd
set(command_line
        ${TEST_PROGRAM}
        --inputDicomDirectory ${INPUTDICOMDIRECTORY}
        --outputVolume ${OUTPUTVOLUME}
        )

message("Running ${command_line}")

execute_process(COMMAND ${command_line}
        RESULT_VARIABLE TEST_RESULT
        )

if(TEST_RESULT)
    message(FATAL_ERROR "Converting SpaceThicknessDiffDicom into Nrrd failed")
endif()

# verify the thickness of nrrd file
set(command_line
        ${KEYVALUE_COMPARE_PROGRAM}
        --f ${OUTPUTVOLUME}
        --tag thicknesses
        --section 2
        --subsection 0
        --v 5
        --numtype 1
        )
message("Running ${command_line}")
execute_process(COMMAND
        ${command_line}
        RESULT_VARIABLE TEST_RESULT
        )
if(TEST_RESULT)
    message(FATAL_ERROR
            "Failed: thickness from nrrd file is not correct")
endif()

# verify the zSpace of nrrd file
set(command_line
        ${KEYVALUE_COMPARE_PROGRAM}
        --f ${OUTPUTVOLUME}
        --tag "space directions"
        --section 2
        --subsection 2
        --v 6
        --numtype 1
        )
message("Running ${command_line}")
execute_process(COMMAND
        ${command_line}
        RESULT_VARIABLE TEST_RESULT
        )
if(TEST_RESULT)
    message(FATAL_ERROR
            "Failed: zSpace from nrrd file is not correct")
endif()

message("Passed")