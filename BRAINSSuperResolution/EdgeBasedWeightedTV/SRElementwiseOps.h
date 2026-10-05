/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef SRElementwiseOps_h_
#define SRElementwiseOps_h_

#include "SRTypes.h"

#include <algorithm>
#include <cblas.h>

// Special override for CVImageType
inline PrecisionType *
GetFirstPointer(CVImageType::Pointer in)
{
  CVType * firstCovariantVectorX = in->GetBufferPointer();
  return firstCovariantVectorX->GetDataPointer();
}

inline PrecisionType *
GetFirstPointer(HalfHermetianImageType::Pointer in)
{
  return reinterpret_cast<PrecisionType *>(in->GetBufferPointer());
}


template <typename ImagePointerType>
PrecisionType *
GetFirstPointer(ImagePointerType in)
{
  return in->GetBufferPointer();
}


// Implement out = c*(a*x + y), y is output, and is corrupted
template <typename ImagePointerType>
void
AddAllElements(ImagePointerType &  OutImg,
               const PrecisionType aScaler,
               ImagePointerType &  xImg,
               ImagePointerType &  yImg,
               const PrecisionType cScaler = 1.0F)
{
  const size_t    N = xImg->GetLargestPossibleRegion().GetNumberOfPixels() * xImg->GetNumberOfComponentsPerPixel();
  PrecisionType * x = GetFirstPointer(xImg);
  PrecisionType * y = GetFirstPointer(yImg);
  cblas_saxpy(N, aScaler, x, 1, y, 1);
  if (cScaler != 1.0F)
  {
    cblas_sscal(N, cScaler, y, 1);
  }
  if (OutImg.GetPointer() != yImg.GetPointer())
  {
    ImagePointerType temp = OutImg;
    OutImg = yImg;
    yImg = temp; // Sanity Check to induce failures if variable is needed in future.
    // yImg has been corrupted by processing here.
  }
}

// Implement out = x*y, y, and y is corrupted
template <typename ImagePointerType>
void
MultiplyVectors(ImagePointerType & OutImg, ImagePointerType & xImg, ImagePointerType & yImg)
{
  const size_t    N = xImg->GetLargestPossibleRegion().GetNumberOfPixels() * xImg->GetNumberOfComponentsPerPixel();
  PrecisionType * x_Start = GetFirstPointer(xImg);
  const PrecisionType * x_End = x_Start + N;
  const PrecisionType * y = GetFirstPointer(yImg);
  PrecisionType *       o = GetFirstPointer(OutImg);
  for (PrecisionType * x = x_Start; x < x_End; ++x)
  {
    (*o) = (*y) * (*x);
    ++o;
    ++y;
  }
}

template <typename ImageTypePointer>
void
Duplicate(ImageTypePointer & Y, ImageTypePointer & YminusL)
{
  const size_t          N = Y->GetLargestPossibleRegion().GetNumberOfPixels() * Y->GetNumberOfComponentsPerPixel();
  const PrecisionType * firstInput = GetFirstPointer(Y);
  const PrecisionType * lastInput = firstInput + N;
  PrecisionType *       firstOutput = GetFirstPointer(YminusL);
  std::copy(firstInput, lastInput, firstOutput);
}


#endif // SRElementwiseOps_h_
