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

#ifndef _gtractDiffusionTensor3D_txx
#define _gtractDiffusionTensor3D_txx

#include "gtractDiffusionTensor3D.h"
#include "itkNumericTraits.h"

namespace itk
{
/**
 *  Compute the Volume Ratio
 */
template <typename TComponent>
typename gtractDiffusionTensor3D<TComponent>::RealValueType
gtractDiffusionTensor3D<TComponent>
::GetVolumeRatio() const
{
  const RealValueType xx = ( *this )[0];
  const RealValueType xy = ( *this )[1];
  const RealValueType xz = ( *this )[2];
  const RealValueType yy = ( *this )[3];
  const RealValueType yz = ( *this )[4];
  const RealValueType zz = ( *this )[5];

  const RealValueType Dav = ( xx + yy + zz ) / 3;

  return 1.0
         - ( ( xx * yy * zz + 2.0
               * ( xy * xz * yz ) - ( zz * xy * xy + yy * xz * xz + xx * yz * yz ) ) / ( Dav * Dav * Dav ) );
}

/**
 *  Compute the Axial Diffusivity
 */
template <typename TComponent>
typename gtractDiffusionTensor3D<TComponent>::RealValueType
gtractDiffusionTensor3D<TComponent>
::GetAxialDiffusivity() const
{
  EigenValuesArrayType eigenValues;

  this->Superclass::ComputeEigenValues( eigenValues );

  // lambda.1
  return vnl_math_abs(eigenValues[0]);
}

/**
 *  Compute the Radial Diffusivity
 */
template <typename TComponent>
typename gtractDiffusionTensor3D<TComponent>::RealValueType
gtractDiffusionTensor3D<TComponent>
::GetRadialDiffusivity() const
{
  EigenValuesArrayType eigenValues;

  this->Superclass::ComputeEigenValues( eigenValues );

  // (lambda.2 + lambda.3)/2.0
  return ( vnl_math_abs(eigenValues[1]) + vnl_math_abs(eigenValues[2]) ) / 2.0;
}

/**
 *  Compute the Lattice Index
 */
template <typename TComponent>
typename gtractDiffusionTensor3D<TComponent>::RealValueType
gtractDiffusionTensor3D<TComponent>
::GetLatticeIndex() const
{
  const RealValueType FA = this->GetFractionalAnisotropy();

  return ( FA + FA * FA ) / 2.0;
}
} // end namespace itk

#endif
