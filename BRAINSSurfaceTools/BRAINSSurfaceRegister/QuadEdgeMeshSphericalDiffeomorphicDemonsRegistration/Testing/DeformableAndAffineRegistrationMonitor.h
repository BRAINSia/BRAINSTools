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

#ifndef __DeformableAndAffineRegistrationMonitor_h
#define __DeformableAndAffineRegistrationMonitor_h

#include "RegistrationMonitor.h"

#include "itkMatrixOffsetTransformBase.h"

#include "vtkMatrix4x4.h"

/** \class DeformableAndAffineRegistrationMonitor
 *  This class provides a VTK visualization pipeline configured for monitoring
 *  the progress of a deformable registration process.
 */
template <class TDeformationFilter, class TPointSet>
class DeformableAndAffineRegistrationMonitor : public RegistrationMonitor
{
public:

  typedef DeformableAndAffineRegistrationMonitor Self;
  typedef RegistrationMonitor                    Superclass;

  typedef TPointSet PointSetType;

  typedef itk::MatrixOffsetTransformBase<double, 3, 3> TransformType;

  DeformableAndAffineRegistrationMonitor();
  virtual ~DeformableAndAffineRegistrationMonitor();

  void ObserveData( const PointSetType * destinationPointSet );

  void ObserveData( const TransformType * transform );

  void ObserveData( const TransformType * transform, const PointSetType * destinationPointSet );

protected:

  virtual void UpdateDataBeforeRendering();

  virtual void PrintOutUpdateMessage();

  typedef enum
    {
    AFFINE,
    DEFORMABLE,
    AFFINEANDDEFORMABLE
    }  RegistrationModeType;

  RegistrationModeType GetRegistrationMode() const;

private:
  void AffineUpdateDataBeforeRendering();

  void DeformableUpdateDataBeforeRendering();

  void AffineAndDeformableUpdateDataBeforeRendering();

  void CopyPointLocationsToFixedSurface();

  void CopyFixedSurfaceMatrix();

  typename PointSetType::ConstPointer   ObservedPointSet;

  vtkSmartPointer<vtkMatrix4x4> Matrix;

  TransformType::ConstPointer ObservedTransform;

  RegistrationModeType RegistrationMode;
};

#include "DeformableAndAffineRegistrationMonitor.txx"

#endif
