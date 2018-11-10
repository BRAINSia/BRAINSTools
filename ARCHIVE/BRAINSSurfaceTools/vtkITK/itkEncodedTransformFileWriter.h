/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEncodedTransformFileWriter.h,v $
  Language:  C++
  Date:      $Date: 2007-08-10 15:43:28 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkEncodedTransformFileWriter_h
#define __itkEncodedTransformFileWriter_h

#include "itkTransformFileWriter.h"
#include "vtkITK.h"

namespace itk
{
class VTK_ITK_EXPORT EncodedTransformFileWriter : public TransformFileWriter
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(EncodedTransformFileWriter);

  /** SmartPointer type alias support */
  using Self = EncodedTransformFileWriter;
  using Pointer = SmartPointer<Self>;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  using Superclass = TransformFileWriter;
  itkTypeMacro(EncodedTransformFileWriter, TransformFileWriter);
  using TransformType = Superclass::TransformType;
  using TransformPointer = Superclass::TransformPointer;

  /** Set/Get the input transform to write */
  void SetInput(const TransformType *transform);

  const TransformType * GetInput() {return *( m_TransformList.begin() ); }

  /** Add a transform to be written */
  void AddTransform(const TransformType *transform);

  /** Write out the transform */
  void Update();

protected:
  EncodedTransformFileWriter();
  virtual ~EncodedTransformFileWriter();
private:
  std::list<const TransformType *> m_TransformList;
};
} // namespace itk

#endif // __itkEncodedTransformFileWriter_h
