project(NetConfiguration)

file(GLOB ProcessDescription_SRC *.cxx)
link_libraries( itksys )
add_library(NetConfiguration ${ProcessDescription_SRC})

install(TARGETS NetConfiguration
                 RUNTIME DESTINATION bin
                 LIBRARY DESTINATION lib
                 ARCHIVE DESTINATION lib/static)

