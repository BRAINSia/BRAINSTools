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

#ifndef _gtractDiffusionTensor3D_hxx
#define _gtractDiffusionTensor3D_hxx

#include "gtractDiffusionTensor3D.h"
#include "itkNumericTraits.h"

namespace itk
{
/**
 *  Compute the Volume Ratio
 */
template <typename TComponent>
typename gtractDiffusionTensor3D<TComponent>::RealValueType
gtractDiffusionTensor3D<TComponent>::GetVolumeRatio() const
{
  const RealValueType xx = (*this)[0];
  const RealValueType xy = (*this)[1];
  const RealValueType xz = (*this)[2];
  const RealValueType yy = (*this)[3];
  const RealValueType yz = (*this)[4];
  const RealValueType zz = (*this)[5];

  const RealValueType Dav = (xx + yy + zz) / 3;

  return 1.0 -
         ((xx * yy * zz + 2.0 * (xy * xz * yz) - (zz * xy * xy + yy * xz * xz + xx * yz * yz)) / (Dav * Dav * Dav));
}

/**
 *  Compute the Axial Diffusivity
 */
template <typename TComponent>
typename gtractDiffusionTensor3D<TComponent>::RealValueType
gtractDiffusionTensor3D<TComponent>::GetAxialDiffusivity() const
{
  EigenValuesArrayType eigenValues;

  this->Superclass::ComputeEigenValues(eigenValues);

  return itk::Math::abs(eigenValues[2]);
}

/**
 *  Compute the Radial Diffusivity
 */
template <typename TComponent>
typename gtractDiffusionTensor3D<TComponent>::RealValueType
gtractDiffusionTensor3D<TComponent>::GetRadialDiffusivity() const
{
  EigenValuesArrayType eigenValues;

  this->Superclass::ComputeEigenValues(eigenValues);

  return (itk::Math::abs(eigenValues[0]) + itk::Math::abs(eigenValues[1])) / 2.0;
}

/**
 *  Compute the Lattice Index
 */
template <typename TComponent>
typename gtractDiffusionTensor3D<TComponent>::RealValueType
gtractDiffusionTensor3D<TComponent>::GetLatticeIndex() const
{
  const RealValueType FA = this->GetFractionalAnisotropy();

  return (FA + FA * FA) / 2.0;
}
} // end namespace itk

#endif
