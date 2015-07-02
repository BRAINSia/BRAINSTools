##- project(AutoWorkup)

#-- Now part of nipype #set(OUTPUT_PATH ${BRAINSTools_CLI_RUNTIME_OUTPUT_DIRECTORY})
#-- Now part of nipype set(OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
#-- Now part of nipype set(OUTPUT_SENTINAL_FILE ${OUTPUT_PATH}/semtools/__init__.py)
#-- Now part of nipype
#-- Now part of nipype message(STATUS "TEST:  Building ${OUTPUT_SENTINAL_FILE}")
#-- Now part of nipype
#-- Now part of nipype add_custom_command(
#-- Now part of nipype    OUTPUT ${OUTPUT_SENTINAL_FILE}
#-- Now part of nipype    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/local_generate_classes.py
#-- Now part of nipype            --python_paths ${BRAINSTools_BINARY_DIR}/../NIPYPE
#-- Now part of nipype            --program_paths ${BRAINSTools_CLI_RUNTIME_OUTPUT_DIRECTORY}
#-- Now part of nipype            --output_path ${OUTPUT_PATH}
#-- Now part of nipype    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/local_generate_classes.py BRAINSFit
#-- Now part of nipype    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#-- Now part of nipype    COMMENT "Generate the semtools python code for nipype")
#-- Now part of nipype
#-- Now part of nipype # Custom target to drive all update commands.
#-- Now part of nipype add_custom_target(semtools DEPENDS ${OUTPUT_SENTINAL_FILE})
#-- Now part of nipype #add_custom_target(${OUTPUT_SENTINAL_FILE} DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/local_generate_classes.py)

#-----------------------------------------------------------------------------
# Output directories.
#
if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
    add_subdirectory(TestSuite)
endif()
