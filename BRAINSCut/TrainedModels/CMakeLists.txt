#
# Install ANN Spherical Domain Image files

#add_subdirectory(2010_SIEMENS)
add_subdirectory(2011_SIEMENS)

file(GLOB AnnProbImages "*.nii.gz")
foreach(IMAGE ${AnnProbImages})
  get_filename_component(image ${IMAGE} NAME)
  configure_file(${IMAGE} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Atlas/ANN_2011_Beta/${image} COPYONLY IMMEDIATE)
endforeach(IMAGE ${AnnProbImages})

install(FILES ${AnnProbImages} DESTINATION bin/Atlas/ANN_2011_Beta)

