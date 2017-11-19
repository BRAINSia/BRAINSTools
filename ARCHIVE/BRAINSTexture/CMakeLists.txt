cmake_minimum_required(VERSION 3.2)
cmake_policy(VERSION 3.2)

project(BRAINSTexture)

find_package(ITK REQUIRED) # NO_MODULE)
include(${ITK_USE_FILE})

add_executable(BRAINSTexture BRAINSTexture.cxx)
target_link_libraries(BRAINSTexture ${ITK_LIBRARIES})

#    1) itkEllipseSpatialObject
#    2) itkSpatialObjectToImageFilter
#    3) itkSimpleContourExtractorImageFilter
