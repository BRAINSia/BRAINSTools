/*
 * Author: Hui Xie    SheenXH@gmail.com
 * date: Sep 20th, 2016
*/

#ifndef BRAINSTOOLS_ITK3DIMAGE2MEXARRAY_H
#define BRAINSTOOLS_ITK3DIMAGE2MEXARRAY_H

#include "mex.h"
#include <string>

extern void itk3DImage2MxArray(const std::string filename, mxArray* plhs[]);

#endif //BRAINSTOOLS_ITK3DIMAGE2MEXARRAY_H
