
#
# Depends on:
#  CMakeParseArguments.cmake from Cmake 2.8.4 or greater
#
include(CMakeParseArguments)

macro(SEMMacroBuildCLI)
  set(options EXECUTABLE_ONLY NO_INSTALL VERBOSE)
  set(oneValueArgs  NAME LOGO_HEADER CLI_XML_FILE CLI_SHARED_LIBRARY_WRAPPER_CXX
                    RUNTIME_OUTPUT_DIRECTORY
                    LIBRARY_OUTPUT_DIRECTORY
                    ARCHIVE_OUTPUT_DIRECTORY
                    INSTALL_RUNTIME_DESTINATION
                    INSTALL_LIBRARY_DESTINATION
                    INSTALL_ARCHIVE_DESTINATION
  )
  set(multiValueArgs ADDITIONAL_SRCS TARGET_LIBRARIES LINK_DIRECTORIES INCLUDE_DIRECTORIES)
  CMAKE_PARSE_ARGUMENTS(LOCAL_SEM
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )


  message(STATUS "Configuring SEM CLI module: ${LOCAL_SEM_NAME}")
  # --------------------------------------------------------------------------
  # Print information helpful for debugging checks
  # --------------------------------------------------------------------------
  if(LOCAL_SEM_VERBOSE)
    list(APPEND ALL_OPTIONS ${options} ${oneValueArgs} ${multiValueArgs})
    foreach(curr_opt ${ALL_OPTIONS})
      message(STATUS "STATUS: ${curr_opt} = ${LOCAL_SEM_${curr_opt}}")
    endforeach()
  endif()
  if(LOCAL_SEM_INSTALL_UNPARSED_ARGUMENTS)
    message(STATUS "WARNING:  Unparsed arguments given [${LOCAL_SEM_INSTALL_UNPARSED_ARGUMENTS}]")
  endif()
  # --------------------------------------------------------------------------
  # Sanity checks
  # --------------------------------------------------------------------------
  if(NOT DEFINED LOCAL_SEM_NAME)
    message(FATAL_ERROR "error: NAME is mandatory: [${LOCAL_SEM_NAME}]")
  endif()

  if(DEFINED LOCAL_SEM_LOGO_HEADER AND NOT EXISTS ${LOCAL_SEM_LOGO_HEADER})
    message(WARNING "warning: Specified LOGO_HEADER [${LOCAL_SEM_LOGO_HEADER}] doesn't exist")
    set(LOCAL_SEM_LOGO_HEADER)
  endif()

  foreach(v LOCAL_SEM_CLI_SHARED_LIBRARY_WRAPPER_CXX)
    if(NOT EXISTS "${${v}}")
      message(FATAL_ERROR "error: Variable ${v} point to an non-existing file or directory !")
    endif()
  endforeach()

  if(DEFINED LOCAL_SEM_CLI_XML_FILE)
    set(cli_xml_file ${LOCAL_SEM_CLI_XML_FILE})
    if(NOT EXISTS ${cli_xml_file})
      message(STATUS "WARNING: Requested Xml file [${cli_xml_file}] doesn't exist !")
    endif()
  else()
    set(cli_xml_file ${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_SEM_NAME}.xml)
    if(NOT EXISTS ${cli_xml_file})
      set(cli_xml_file ${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_SEM_NAME}.xml)
      if(NOT EXISTS ${cli_xml_file})
        message(STATUS "WARNING: Default SOURECE_DIR Xml file [${cli_xml_file}] doesn't exist !")
        message(STATUS "WARNING: Default BINARY_DIR Xml file [${cli_xml_file}] doesn't exist !")
      endif()
    endif()
  endif()
  if(NOT EXISTS ${cli_xml_file})
      message(FATAL_ERROR "Xml file [${cli_xml_file}] doesn't exist !")
  endif()

  set(CLP ${LOCAL_SEM_NAME})

  # SlicerExecutionModel
  find_package(SlicerExecutionModel REQUIRED GenerateCLP)
  include(${GenerateCLP_USE_FILE})

  set(${CLP}_SOURCE ${CLP}.cxx ${LOCAL_SEM_ADDITIONAL_SRCS})
  generateclp(${CLP}_SOURCE ${cli_xml_file} ${LOCAL_SEM_LOGO_HEADER})

  if(DEFINED LOCAL_SEM_LINK_DIRECTORIES)
    link_directories(${LOCAL_SEM_LINK_DIRECTORIES})
  endif()

  if(DEFINED LOCAL_SEM_INCLUDE_DIRECTORIES)
    include_directories(${LOCAL_SEM_INCLUDE_DIRECTORIES})
  endif()

  set(cli_targets)

  if(NOT LOCAL_SEM_EXECUTABLE_ONLY)

    add_library(${CLP}Lib SHARED ${${CLP}_SOURCE})
    set_target_properties(${CLP}Lib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
    if(DEFINED LOCAL_SEM_TARGET_LIBRARIES)
      target_link_libraries(${CLP}Lib ${LOCAL_SEM_TARGET_LIBRARIES})
    endif()

    add_executable(${CLP} ${LOCAL_SEM_CLI_SHARED_LIBRARY_WRAPPER_CXX})
    target_link_libraries(${CLP} ${CLP}Lib)

    set(cli_targets ${CLP} ${CLP}Lib)

  else()

    add_executable(${CLP} ${${CLP}_SOURCE})
    if(DEFINED LOCAL_SEM_TARGET_LIBRARIES)
      target_link_libraries(${CLP} ${LOCAL_SEM_TARGET_LIBRARIES})
    endif()

    set(cli_targets ${CLP})

  endif()

  # Set labels associated with the target.
  set_target_properties(${cli_targets} PROPERTIES LABELS ${CLP})

  if(NOT DEFINED LOCAL_SEM_RUNTIME_OUTPUT_DIRECTORY)
    set(LOCAL_SEM_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    message(STATUS "WARNING:  Setting default RUNTIME_OUTPUT_DIRECTORY to ${LOCAL_SEM_RUNTIME_OUTPUT_DIRECTORY}")
  endif()
  if(NOT DEFINED LOCAL_SEM_LIBRARY_OUTPUT_DIRECTORY)
    set(LOCAL_SEM_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    message(STATUS "WARNING:  Setting default LIBRARY_OUTPUT_DIRECTORY to ${LOCAL_SEM_LIBRARY_OUTPUT_DIRECTORY}")
  endif()
  if(NOT DEFINED LOCAL_SEM_ARCHIVE_OUTPUT_DIRECTORY)
    set(LOCAL_SEM_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    message(STATUS "WARNING:  Setting default ARCHIVE_OUTPUT_DIRECTORY to ${LOCAL_SEM_ARCHIVE_OUTPUT_DIRECTORY}")
  endif()
  if(NOT DEFINED LOCAL_SEM_INSTALL_RUNTIME_DESTINATION)
    set(LOCAL_SEM_INSTALL_RUNTIME_DESTINATION ${CMAKE_BINARY_DIR}/bin)
    message(STATUS "WARNING:  Setting default INSTALL_RUNTIME_DESTINATION to ${LOCAL_SEM_INSTALL_RUNTIME_DESTINATION}")
  endif()
  if(NOT DEFINED LOCAL_SEM_INSTALL_LIBRARY_DESTINATION)
    set(LOCAL_SEM_INSTALL_LIBRARY_DESTINATION ${CMAKE_BINARY_DIR}/lib)
    message(STATUS "WARNING:  Setting default INSTALL_LIBRARY_DESTINATION to ${LOCAL_SEM_INSTALL_LIBRARY_DESTINATION}")
  endif()
  if(NOT DEFINED LOCAL_SEM_INSTALL_ARCHIVE_DESTINATION)
    set(LOCAL_SEM_INSTALL_ARCHIVE_DESTINATION ${CMAKE_BINARY_DIR}/lib)
    message(STATUS "WARNING:  Setting default INSTALL_ARCHIVE_DESTINATION to ${LOCAL_SEM_INSTALL_ARCHIVE_DESTINATION}")
  endif()

  set_target_properties(${cli_targets} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${LOCAL_SEM_RUNTIME_OUTPUT_DIRECTORY}"
    LIBRARY_OUTPUT_DIRECTORY "${LOCAL_SEM_RUNTIME_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY "${LOCAL_SEM_RUNTIME_OUTPUT_DIRECTORY}"
    )

  if(NOT LOCAL_SEM_NO_INSTALL)
    # Install each target in the production area (where it would appear in an installation)
    # and install each target in the developer area (for running from a build)
    install(TARGETS ${cli_targets}
      RUNTIME DESTINATION ${LOCAL_SEM_INSTALL_RUNTIME_DESTINATION} COMPONENT RuntimeLibraries
      LIBRARY DESTINATION ${LOCAL_SEM_INSTALL_LIBRARY_DESTINATION} COMPONENT RuntimeLibraries
      ARCHIVE DESTINATION ${LOCAL_SEM_INSTALL_ARCHIVE_DESTINATION} COMPONENT RuntimeLibraries
      )
  endif()


endmacro()
