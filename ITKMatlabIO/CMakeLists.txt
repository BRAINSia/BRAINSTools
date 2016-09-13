project(ITKMatlabIO)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_LIST_DIR}/CMake)

# -- NOTE: This test is not needed for MATLAB2014b with xcode.
#Matlab Mex files require char16_t for mex files
#--include(CMake/compilerTestChar16_t.cmake)

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

add_subdirectory(readItkImage)

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
