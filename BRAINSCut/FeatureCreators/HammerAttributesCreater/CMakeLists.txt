project(HammerAttributeCreator)

cmake_minimum_required(VERSION 2.8)

find_package(ITK)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
else(ITK_FOUND)
  message(FATAL_ERROR "ITK not found.  Please set ITK_DIR.")
endif(ITK_FOUND)

set(CLP HammerAttributeCreator)
set(${CLP}_SOURCE ${CLP}.cxx)
generateclp(${CLP}_SOURCE ${CLP}.xml ${BRAINSLogo})

# Build it
add_executable(${CLP} ${${CLP}_SOURCE})
#slicer3_set_plugins_output_path(${CLP})
target_link_libraries(${CLP} ${ITK_LIBRARIES})

