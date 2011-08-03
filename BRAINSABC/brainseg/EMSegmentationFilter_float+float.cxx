#include "EMSegmentationFilter.txx"

#include "itkImage.h"

typedef itk::Image<float, 3> FloatImageType;

template class EMSegmentationFilter<FloatImageType, FloatImageType>;
