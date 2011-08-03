#ifndef __filterFloatImages_h
#define -_filterFloatImages_h

#include "itkImage.h"

#include <vector>

#include <string>

void filterFloatImages(std::vector<itk::Image<float, 3>::Pointer> & images, std::string & method, unsigned int iters,
                       double dt);

#endif
