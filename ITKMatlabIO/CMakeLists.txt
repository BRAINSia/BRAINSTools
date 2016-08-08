project(ITKMatlabIO)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMake)

# -- NOTE: This test is not needed for MATLAB2014b with xcode.
#Matlab Mex files require char16_t for mex files
#--include(CMake/compilerTestChar16_t.cmake)

option(${PRIMARY_PROJECT_NAME}_INSTALL_DEVELOPMENT "Install development support include and libraries for external packages." OFF)
mark_as_advanced(${PRIMARY_PROJECT_NAME}_INSTALL_DEVELOPMENT)

# Some test are failing due to inadequate test construction, but
# the code seems to do the correct thing on real data.
# This is also for tests that are huge, and can not be running
# a long time.
option(ENABLE_EXTENDED_TESTING "Enable tests that are long running, or where the test itself is in error." OFF)
mark_as_advanced(ENABLE_EXTENDED_TESTING)

#-----------------------------------------------------------------------
set(TestData_DIR ${CMAKE_CURRENT_SOURCE_DIR}/TestData)

#-----------------------------------------------------------------------------
# Define list of module names
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Add module sub-directory if USE_<MODULENAME> is both defined and true
#-----------------------------------------------------------------------------
add_subdirectory(mexFiles)

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
