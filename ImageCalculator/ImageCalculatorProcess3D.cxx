/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    ImageCalculatorProcess3D.cxx
  Language:  C++

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "ImageCalculatorTemplates.h"

// This program calls the 3d function to process the algorithm.
void ImageCalculatorProcess3D(const std::string & InType, MetaCommand & command)
{
  ImageCalculatorProcessND<3>(InType, command);
}
