project(ProcessDescription)

file(GLOB ProcessDescription_SRC *.cxx)
link_libraries( itksys )
add_library(ProcessDescription ${ProcessDescription_SRC})

install(TARGETS ProcessDescription
                 RUNTIME DESTINATION bin
                 LIBRARY DESTINATION lib
                 ARCHIVE DESTINATION lib/static)

