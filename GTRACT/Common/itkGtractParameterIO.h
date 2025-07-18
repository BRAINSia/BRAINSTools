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
/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkGtractParameterIO_h
#define __itkGtractParameterIO_h

#include <iostream>

#include "itkObject.h"
#include "itkImage.h"
#include "itkRGBAPixel.h"

#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaImageIO.h"
#include "gtractCommonWin32.h"

#include "algo.h"

namespace itk
{
/** \class GtractImageIO
 * \brief Convience functions for DTI parameter files. These were required
 * for building on Windows Visual studio because the object sizes
 * from the templated code.
 */
class GTRACT_COMMON_EXPORT GtractParameterIO : public itk::Object
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(GtractParameterIO);

  using Self = GtractParameterIO;
  using Superclass = itk::Object;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;
  itkOverrideGetNameOfClassMacro(gtractParameterIO);
  itkNewMacro(Self);

  /*** Get/Set the Images for I/O Routines ***/
  itkGetMacro(Directions, TMatrix);
  itkGetMacro(Bvalue, float);
  itkGetMacro(NumberOfDirections, int);
  itkGetMacro(FileName, std::string);
  itkSetMacro(FileName, std::string);

  /*** Additional API - Functions ***/
  void
  Update();

protected:
  /** Constructor and Destructor */
  GtractParameterIO();
  ~GtractParameterIO() override = default;

private:
  TMatrix     m_Directions;
  float       m_Bvalue{ 0.0 };
  int         m_NumberOfDirections{ 0 };
  std::string m_FileName;
};
} // namespace itk
#endif
