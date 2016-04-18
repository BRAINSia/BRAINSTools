#ifndef BRAINSTypes_H__
#define BRAINSTypes_H__

#include "itkImageMaskSpatialObject.h"
#include <map>
#include <string>

typedef itk::SpatialObject<3>      SpatialObjectType;
typedef SpatialObjectType::Pointer ImageMaskPointer;

typedef itk::Point<double, 3>                    LandmarkPointType;
typedef std::map<std::string, LandmarkPointType> LandmarksMapType;
typedef std::map<std::string, float>             LandmarkWeightMapType;
typedef std::map<std::string, double>            LandmarksWeightMapType;

#endif //BRAINSTypes_H__
