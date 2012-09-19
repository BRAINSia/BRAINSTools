##- project(BRAINSCutConfiguration)

file(GLOB NetConfiguration_SRCS *.cxx)
link_libraries( itksys )

##
## Build NetConfigurationCOMMONLIB library
##
add_library(NetConfigurationCOMMONLIB STATIC ${NetConfiguration_SRCS})
target_link_libraries(NetConfigurationCOMMONLIB BRAINSCommonLib ${ITK_LIBRARIES})
