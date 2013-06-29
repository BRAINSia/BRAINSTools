##- project(AutoWorkup)

#set(OUTPUT_PATH ${BRAINSTools_CLI_RUNTIME_OUTPUT_DIRECTORY})
set(OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(OUTPUT_SENTINAL_FILE ${OUTPUT_PATH}/SEMTools/__init__.py)

message(STATUS "TEST:  Building ${OUTPUT_SENTINAL_FILE}")

add_custom_command(
   OUTPUT ${OUTPUT_SENTINAL_FILE}
   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/local_generate_classes.py
           --python_paths ${BRAINSTools_CLI_RUNTIME_OUTPUT_DIRECTORY}/../NIPYPE
           --program_paths ${BRAINSTools_CLI_RUNTIME_OUTPUT_DIRECTORY}
           --output_path ${OUTPUT_PATH}
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/local_generate_classes.py
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   COMMENT "Generate the SEMTools python code for nipype")

# Custom target to drive all update commands.
add_custom_target(SEMTools DEPENDS ${OUTPUT_SENTINAL_FILE})
#add_custom_target(${OUTPUT_SENTINAL_FILE} DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/local_generate_classes.py)

#-----------------------------------------------------------------------------
# Output directories.
#
if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
    add_subdirectory(TestSuite)
endif()
