#
# Building external projects via ExternalProject requires passing along some important CMkae variables into the autoconf configuration.
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

macro(AutoConf_FLAGS varname lang additional_flags)
  string(TOUPPER CONF ${CMAKE_BUILD_TYPE})
  set(${varname} "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${CONF}} ${additional_flags}")
endmacro()