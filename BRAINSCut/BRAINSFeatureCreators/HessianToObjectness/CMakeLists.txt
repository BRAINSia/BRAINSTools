##- project(CannyEdge)

## Build all the programs
##
add_executable(HessianToObjectness HessianToObjectness.cxx)
target_link_libraries(HessianToObjectness BRAINSCommonLib)
