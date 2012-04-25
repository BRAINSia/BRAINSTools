/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTestingMacros.h,v $
  Language:  C++
  Date:      $Date: 2009-05-09 17:40:20 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __AffineRegistrationMonitor_h
#define __AffineRegistrationMonitor_h

#include "RegistrationMonitor.h"
#include "itkMatrixOffsetTransformBase.h"

#include "vtkMatrix4x4.h"

/** \class AffineRegistrationMonitor
 *  This class provides a VTK visualization pipeline configured for monitoring
 *  the progress of a registration process.
 */
class AffineRegistrationMonitor : public RegistrationMonitor
{
public:

  typedef AffineRegistrationMonitor Self;
  typedef RegistrationMonitor       Superclass;

  typedef itk::MatrixOffsetTransformBase<double, 3, 3> TransformType;

  AffineRegistrationMonitor();
  virtual ~AffineRegistrationMonitor();

  void ObserveData( const TransformType * transform );

protected:

  virtual void UpdateDataBeforeRendering();

private:

  vtkSmartPointer<vtkMatrix4x4> Matrix;

  TransformType::ConstPointer ObservedTransform;
};

#endif
