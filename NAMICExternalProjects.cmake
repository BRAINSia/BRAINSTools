
include(${CMAKE_CURRENT_LIST_DIR}/Common.cmake)

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
  ${${PROJECT_NAME}_SOURCE_DIR}/CMake
  ${${PROJECT_NAME}_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )

#-----------------------------------------------------------------------------
find_package(ITK REQUIRED)
if(Slicer_BUILD_${PROJECT_NAME})
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1) //Incorporate with Slicer Nicely
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${GenerateCLP_USE_FILE})
include(${SlicerExecutionModel_USE_FILE})
include(${SlicerExecutionModel_CMAKE_DIR}/SEMMacroBuildCLI.cmake)

#-----------------------------------------------------------------------------
enable_testing()
include(CTest)

#-----------------------------------------------------------------------
# Setup locations to find externally maintained test data.
#-----------------------------------------------------------------------
include(${PROJECT_NAME}ExternalData)

add_subdirectory(src)

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
