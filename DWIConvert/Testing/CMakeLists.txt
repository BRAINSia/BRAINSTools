#
# Test CLI executables
#
include_regular_expression("^.*$")

# SlicerExecutionModel
find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${GenerateCLP_USE_FILE})

#include_directories(${ModuleDescriptionParser_SOURCE_DIR})
#include_directories (${TCLAP_SOURCE_DIR}/include)

set (CLP DWIConvert)
add_executable(${CLP}Test ${CLP}Test.cxx)
target_include_directories(${CLP}Test PRIVATE
  ${DWIConvert_SOURCE_DIR}
  ${DWIConvert_BINARY_DIR}
)
add_dependencies(${CLP}Test ${CLP})

target_link_libraries(${CLP}Test PRIVATE ${CLP}Lib ${LIBICONV_LIBRARIES})

set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})
set_target_properties(${CLP}Test PROPERTIES FOLDER ${MODULE_FOLDER})

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------

set(TEST_DATA "${Slicer_SOURCE_DIR}/Testing/Data/Input")
set(TEMP "${Slicer_BINARY_DIR}/Testing/Temporary")

# DWIConvert tests
add_test(
    NAME ${CLP}Test
    COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
    ModuleEntryPoint
    --useBMatrixGradientDirections
    --inputDicomDirectory ${TEST_DATA}/DWIDicom
    --outputDirectory ${TEMP}
    --outputVolume DWIDicom.nhdr
  )
set_property(TEST ${CLP}Test PROPERTY LABELS ${CLP})
