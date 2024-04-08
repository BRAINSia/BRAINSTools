# Dicom file has different zSpace and Thickness data test
# first, convert dicom data into Nrrd
set(command_line
        ${TEST_PROGRAM}
        --inputDicomDirectory ${INPUTDICOMDIRECTORY}
        --outputVolume ${OUTPUTVOLUME}
        )

message("Running ${command_line}")

execute_process(COMMAND ${command_line}
        RESULT_VARIABLE TEST_RESULT1
        )

if(TEST_RESULT1)
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
        RESULT_VARIABLE TEST_RESULT2
        )
if(TEST_RESULT2)
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
        RESULT_VARIABLE TEST_RESULT3
        )
if(TEST_RESULT3)
    message(FATAL_ERROR
            "Failed: zSpace from nrrd file is not correct")
endif()

if (TEST_RESULT1 OR TEST_RESULT2 OR TEST_RESULT3)
    message("Failed: checking the case that zSpace and thickness have different value.")
else()
    message("Passed")
endif()
