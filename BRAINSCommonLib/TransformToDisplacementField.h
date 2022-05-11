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
#ifndef __TransformToDisplacementField_h
#define __TransformToDisplacementField_h

#include "itkIO.h"
#include "CrossOverAffineSystem.h"

#include <itkTransformToDisplacementFieldFilter.h>

/**
 * Go from any subclass of Transform, to the corresponding deformation field
 */
template <typename DisplacementFieldPointerType, typename TransformPointerType>
DisplacementFieldPointerType
TransformToDisplacementField(itk::ImageBase<DisplacementFieldPointerType::ObjectType::ImageDimension> * templateImage,
                             TransformPointerType                                                       xfrm)
{
  using OutputType = typename DisplacementFieldPointerType::ObjectType;
  using TodefType = typename itk::TransformToDisplacementFieldFilter<OutputType, double>;
  typename TodefType::Pointer todef(TodefType::New());
  todef->SetUseReferenceImage(true);
  todef->SetReferenceImage(templateImage);
  todef->SetTransform(xfrm);
  try
  {
    todef->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    throw; // pass the buck up.
  }
  return todef->GetOutput();
}

#endif // TransformToDisplacementField_h
