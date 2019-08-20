get_filename_component(_ITKExternalData_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

#Need a long timeout for DWI huge data sets
set(ExternalData_TIMEOUT_INACTIVITY 900)
set(ExternalData_TIMEOUT_ABSOLUTE   1360)

if(NOT COMMAND ExternalData_add_test)
  if(EXISTS ${CMAKE_ROOT}/Modules/ExternalData.cmake)
    include(ExternalData)
  else()
    include(${_ITKExternalData_DIR}/ExternalData.cmake)
  endif()
endif()


## The user can specify an environmental variable for shared locations of
## of object files to prevent continous downloading of common objects.
if(NOT ExternalData_OBJECT_STORES)
  # Use ExternalData_OBJECT_STORES from environment as default.
  set(ExternalData_OBJECT_STORES_DEFAULT "")
  if(DEFINED "ENV{ExternalData_OBJECT_STORES}")
    file(TO_CMAKE_PATH "$ENV{ExternalData_OBJECT_STORES}" ExternalData_OBJECT_STORES_DEFAULT)
  endif()
endif()

set(ExternalData_OBJECT_STORES "${ExternalData_OBJECT_STORES_DEFAULT}" CACHE STRING
  "Semicolon-separated list of local directories holding data objects in the layout %(algo)/%(hash).")
mark_as_advanced(ExternalData_OBJECT_STORES)
if(NOT ExternalData_OBJECT_STORES)
  set(ExternalData_OBJECT_STORES "${CMAKE_BINARY_DIR}/ExternalData/Objects")
  file(MAKE_DIRECTORY "${ExternalData_OBJECT_STORES}")
endif()
list(APPEND ExternalData_OBJECT_STORES
  # Local data store populated by the ITK pre-commit hook
  "${CMAKE_SOURCE_DIR}/.ExternalData"
  )

set(ExternalData_BINARY_ROOT ${CMAKE_BINARY_DIR}/ExternalData)

set(ExternalData_URL_TEMPLATES "" CACHE STRING
  "Additional URL templates for the ExternalData CMake script to look for testing data. E.g.
file:///var/bigharddrive/%(algo)/%(hash)")
mark_as_advanced(ExternalData_URL_TEMPLATES)
list(APPEND ExternalData_URL_TEMPLATES
  # Local data store populated by the ITK pre-commit hook
  "file:///${${PROJECT_NAME}_SOURCE_DIR}/.ExternalData/%(algo)/%(hash)"
  # Data published by Iowa Psychiatry web interface
  ## The primary home for data
  "http://slicer.kitware.com/midas3/api/rest?method=midas.bitstream.download&checksum=%(hash)"
  # Data published by MIDAS
  "http://midas3.kitware.com/midas/api/rest?method=midas.bitstream.download&checksum=%(hash)&algorithm=%(algo)"

  # Data published by developers using git-gerrit-push.
  "http://www.itk.org/files/ExternalData/%(algo)/%(hash)"
  )

# Tell ExternalData commands to transform raw files to content links.
# INFO: Condition this feature on presence of our pre-commit hook.
set(ExternalData_LINK_CONTENT SHA512)

# Match series of the form <base>.<ext>, <base>.<n>.<ext> such that <base> may
# end in a (test) number that is not part of any series numbering.
set(ExternalData_SERIES_PARSE "()(\\.[^./]*)$")
set(ExternalData_SERIES_MATCH "(\\.[0-9]+)?")
