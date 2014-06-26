##- project(BRAINSCutConfiguration)

file(GLOB NetConfiguration_SRCS *.cxx)
link_libraries( itksys )

##
## Build NetConfigurationCOMMONLIB library
##
add_library(NetConfigurationCOMMONLIB STATIC ${NetConfiguration_SRCS})
target_link_libraries(NetConfigurationCOMMONLIB BRAINSCommonLib ${ITK_LIBRARIES})

#
# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
#
# XXX When minimum CMake version will be changed to version >=2.8.9, the following code
# could be used instead:
# set_target_properties(NetConfigurationCOMMONLIB PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set_target_properties(NetConfigurationCOMMONLIB PROPERTIES COMPILE_FLAGS "-fPIC")
endif()
