/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEncodedTransformFileReader.h,v $
  Language:  C++
  Date:      $Date: 2007-08-09 15:15:39 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkEncodedTransformFileReader_h
#define __itkEncodedTransformFileReader_h

#ifdef ITK_USE_TRANSFORM_IO_FACTORIES
#  include "itkEncodedTransformFileReaderWithFactory.h"
#else

#  include "itkTransformFileReader.h"
#  include "vtkITK.h"

namespace itk
{
class VTK_ITK_EXPORT EncodedTransformFileReader : public TransformFileReader
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(EncodedTransformFileReader);

  /** SmartPointer type alias support */
  using Self = EncodedTransformFileReader;
  using Pointer = SmartPointer<Self>;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  using Superclass = TransformFileReader;
  itkTypeMacro(EncodedTransformFileReader, TransformFileReader);
  using TransformType = Superclass::TransformType;
  using ParametersType = Superclass::ParametersType;
  using TransformPointer = Superclass::TransformPointer;
  using TransformListType = Superclass::TransformListType;

  /** Read in the transform */
  void
  Update();

protected:
  EncodedTransformFileReader();
  virtual ~EncodedTransformFileReader();
};
} // namespace itk

#endif

#endif // __itkEncodedTransformFileReader_h
