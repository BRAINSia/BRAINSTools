#
# Install Trained ANN files

file(GLOB AnnWeightFiles "*.txtANN000000001")
foreach(WEIGHT ${AnnWeightFiles})
  get_filename_component(weight ${WEIGHT} NAME)
  configure_file(${WEIGHT} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Atlas/ANN_2010_Beta/${weight} COPYONLY IMMEDIATE)
endforeach(WEIGHT ${AnnWeightFiles})

install(FILES ${AnnWeightFiles} DESTINATION bin/Atlas/ANN_2010_Beta)

#
# Install ANN Probability Maps files

file(GLOB AnnProbImages "*.nii.gz")
foreach(IMAGE ${AnnProbImages})
  get_filename_component(image ${IMAGE} NAME)
  configure_file(${IMAGE} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Atlas/ANN_2010_Beta/${image} COPYONLY IMMEDIATE)
endforeach(IMAGE ${AnnProbImages})

install(FILES ${AnnProbImages} DESTINATION bin/Atlas/ANN/ANN_2010_Beta)

