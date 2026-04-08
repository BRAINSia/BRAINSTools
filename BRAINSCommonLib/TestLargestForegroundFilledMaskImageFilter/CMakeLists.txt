##The start of a test program for this filter.
## HACK-- This still needs to be finalized, and tests written.

set(TestName TestLargestForegroundFilledMaskImageFilter)

add_executable(${TestName} ${TestName}.cxx)
target_link_libraries(${TestName} PRIVATE BRAINSCommonLib)
set_target_properties(${TestName} PROPERTIES FOLDER ${MODULE_FOLDER})
