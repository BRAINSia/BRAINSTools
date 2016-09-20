/*
 * Author: Hui Xie    SheenXH@gmail.com
 * date: Sep 20th, 2016
*/

#ifndef ITK_IMAGE_2_MEX_ARRAY_H
#define ITK_IMAGE_2_MEX_ARRAY_H

#include "mex.h"
#include <string>

extern void itkImage2MxArray(const std::string filename, mxArray* plhs[]);

#endif //BRAINSTOOLS_CONVERTIMAGE_H
