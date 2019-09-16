#ifndef BRAINSTypes_H__
#define BRAINSTypes_H__

#include "itkImageMaskSpatialObject.h"
#include <map>
#include <string>

using SpatialObjectType = itk::SpatialObject<3>;
using ImageMaskPointer = SpatialObjectType::Pointer;

using LandmarkPointType = itk::Point<double, 3>;
using LandmarksMapType = std::map<std::string, LandmarkPointType>;
// using LandmarksWeightMapType = std::map<std::string, float>;
using LandmarksWeightMapType = std::map<std::string, double>;

#endif // BRAINSTypes_H__
