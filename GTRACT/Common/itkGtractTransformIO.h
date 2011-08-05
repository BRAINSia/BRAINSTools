#error "This file should not be used"
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

#ifndef __itkGtractTransformIO_h
#define __itkGtractTransformIO_h

#include "itkObject.h"
#include "itkBSplineDeformableTransform.h"
#include "itkIdentityTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkAffineTransform.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"
#include "gtractCommonWin32.h"

#include "GenericTransformImage.h"

namespace itk
{
/** \class GtractImageIO
 * \brief Convience functions for transform I/O. These were required
 * for building on Windows Visual studio because the object sizes
 * from the templated code.
 */
class GTRACT_COMMON_EXPORT GtractTransformIO : public itk::Object
{
public:
  typedef GtractTransformIO             Self;
  typedef itk::Object                   Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  itkTypeMacro(GtractTransformIO, itk::Object);
  itkNewMacro(Self);

  typedef double BSplineCoordinateRepType;
  typedef itk::BSplineDeformableTransform<
      BSplineCoordinateRepType, 3, 3>     BSplineTransformType;

  typedef itk::VersorRigid3DTransform<double>                                    RigidTransformType;
  typedef itk::AffineTransform<double, 3>                                        AffineTransformType;
  typedef itk::IdentityTransform<double, 3>                                      IdentityTransformType;
  typedef itk::ThinPlateR2LogRSplineKernelTransform<BSplineCoordinateRepType, 3> ThinPlateSplineTransformType;

  /*** Get/Set the Images for I/O Routines ***/
  itkSetObjectMacro(RigidTransform, RigidTransformType);
  itkSetObjectMacro(AffineTransform, AffineTransformType);
  itkSetObjectMacro(BSplineTransform, BSplineTransformType);
  itkSetObjectMacro(InverseBSplineTransform, ThinPlateSplineTransformType);

  itkGetObjectMacro(RigidTransform, RigidTransformType);
  itkGetObjectMacro(AffineTransform, AffineTransformType);
  itkGetObjectMacro(BSplineTransform, BSplineTransformType);
  itkGetObjectMacro(InverseBSplineTransform, ThinPlateSplineTransformType);

  /*** Additional API - Functions ***/
  void LoadTransform();

  void SaveTransform( int type );

  void SetFileName(char *);

  void SetFileName(std::string);

protected:

  /** Constructor and Destructor */
  GtractTransformIO();
  ~GtractTransformIO()
  {
  }

private:
  GtractTransformIO( const Self & );   // purposely not implemented
  void operator=( const Self & );      // purposely not implemented

  RigidTransformType::Pointer           m_RigidTransform;
  AffineTransformType::Pointer          m_AffineTransform;
  BSplineTransformType::Pointer         m_BSplineTransform;
  ThinPlateSplineTransformType::Pointer m_InverseBSplineTransform;

  std::string m_FileName;
};
}

#endif
