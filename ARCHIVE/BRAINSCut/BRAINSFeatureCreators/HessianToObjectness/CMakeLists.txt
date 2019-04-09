##- project(CannyEdge)

set(MODULE_FOLDER "Module-HessianToObjectness")

## Build all the programs
##
add_executable(HessianToObjectness HessianToObjectness.cxx)
target_link_libraries(HessianToObjectness BRAINSCommonLib)
set_target_properties(HessianToObjectness PROPERTIES FOLDER ${MODULE_FOLDER})
