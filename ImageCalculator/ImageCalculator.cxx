/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    ImageCalculator.cxx
  Language:  C++

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "ImageCalculatorUtils.h"
#ifdef ITKV3_COMPATIBILITY
#include "itkAnalyzeImageIOFactory.h"
#endif
int main(int argc, char *argv[])
{
#ifdef ITKV3_COMPATIBILITY
  itk::ObjectFactoryBase::RegisterFactory( itk::AnalyzeImageIOFactory::New() );
#endif
  return PrimaryImageCalculatorRoutine(argc, argv);
}

// End of main()
