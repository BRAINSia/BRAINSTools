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
  typedef GtractParameterIO             Self;
  typedef itk::Object                   Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  itkTypeMacro(gtractParameterIO, itk::Object);
  itkNewMacro(Self);

  /*** Get/Set the Images for I/O Routines ***/
  itkGetMacro(Directions, TMatrix);
  itkGetMacro(Bvalue, float);
  itkGetMacro( NumberOfDirections,            int );
  itkGetMacro(FileName, std::string);
  itkSetMacro(FileName, std::string);

  /*** Additional API - Functions ***/
  void Update();

protected:
  /** Constructor and Destructor */
  GtractParameterIO();
  ~GtractParameterIO()
  {
  }

private:
  GtractParameterIO( const Self & );     // purposely not implemented
  void operator=( const Self & );        // purposely not implemented

  TMatrix     m_Directions;
  float       m_Bvalue;
  int         m_NumberOfDirections;
  std::string m_FileName;
};
}
#endif
