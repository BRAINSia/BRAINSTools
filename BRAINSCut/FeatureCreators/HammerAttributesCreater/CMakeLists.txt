project(HammerAttributeCreator)

set(CLP HammerAttributeCreator)
set(${CLP}_SOURCE ${CLP}.cxx)
generateclp(${CLP}_SOURCE ${CLP}.xml ${BRAINSLogo})

# Build it
add_executable(${CLP} ${${CLP}_SOURCE})
#slicer3_set_plugins_output_path(${CLP})
target_link_libraries(${CLP} ${ITK_LIBRARIES})

