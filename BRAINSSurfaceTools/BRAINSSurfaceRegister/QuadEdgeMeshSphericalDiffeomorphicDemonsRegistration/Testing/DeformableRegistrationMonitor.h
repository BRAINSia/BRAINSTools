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

#ifndef __DeformableRegistrationMonitor_h
#define __DeformableRegistrationMonitor_h

#include "RegistrationMonitor.h"

/** \class DeformableRegistrationMonitor
 *  This class provides a VTK visualization pipeline configured for monitoring
 *  the progress of a deformable registration process.
 */
template <class TDeformationFilter, class TPointSet>
class DeformableRegistrationMonitor : public RegistrationMonitor
{
public:

  typedef DeformableRegistrationMonitor Self;
  typedef RegistrationMonitor           Superclass;

  typedef TPointSet PointSetType;

  DeformableRegistrationMonitor();
  virtual ~DeformableRegistrationMonitor();

  void ObserveData( const PointSetType * destinationPointSet );

protected:

  virtual void UpdateDataBeforeRendering();

  virtual void PrintOutUpdateMessage();

private:

  typename PointSetType::ConstPointer     ObservedPointSet;
};

#include "DeformableRegistrationMonitor.txx"

#endif
