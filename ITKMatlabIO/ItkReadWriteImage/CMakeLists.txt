#
# Author: Hui Xie    SheenXH@gmail.com
# date: Sep 20th, 2016
#

find_package(ITK)
find_package(Matlab COMPONENTS MX_LIBRARY ENG_LIBRARY REQUIRED)

include_directories(${ITK_INCLUDE_DIRS})
include_directories(${Matlab_INCLUDE_DIRS})

#add_executable(itkReadImage itkReadImage.cpp convertItkImageToMxArray.cpp)
#target_link_libraries(itkReadImage ${ITK_LIBRARIES} ${VTK_LIBRARIES})

#####################
#####################
# Provide a macro to build the mex files from
# within CMake
#####################
include(CMakeParseArguments)
include(${ITK_USE_FILE})

# BuildMex -- arguments
# MEXNAME = root of mex library name
# TARGETDIR = location for the mex library files to be created
# SOURCE = list of source files
# LIBRARIES = libraries needed to link mex library
macro(BuildMex)
    set(oneValueArgs MEXNAME TARGETDIR)
    set(multiValueArgs SOURCE LIBRARIES)
    cmake_parse_arguments(BuildMex "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # message("MEXNAME=${BuildMex_MEXNAME} SOURCE=${BuildMex_SOURCE} LIBRARIES=${BuildMex_LIBRARIES}")
    #--  set_source_files_properties(${BuildMex_SOURCE} COMPILE_DEFINITIONS -DMatlab_MEX_FILE )
    matlab_add_mex(
            NAME ${BuildMex_MEXNAME}
            SRC ${BuildMex_SOURCE}
            #--      [OUTPUT_NAME output_name]
            #--      [DOCUMENTATION file.txt]
            LINK_TO ${BuildMex_LIBRARIES}
            #--      [...]
    )
endmacro(BuildMex)

matlab_add_mex(
        NAME itkReadImage
        SRC itkReadImage.cpp itkImageToMxArray.cpp
        #--      [DOCUMENTATION file.txt]
        LINK_TO ${ITK_LIBRARIES} ${VTK_LIBRARIES}
)
