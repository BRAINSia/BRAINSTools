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

#ifndef __gtractDiffusionTensor3D_h
#define __gtractDiffusionTensor3D_h

// Undefine an eventual DiffusionTensor3D macro
#ifdef DiffusionTensor3D
#  undef DiffusionTensor3D
#endif

#include "itkDiffusionTensor3D.h"

namespace itk
{
/** \class DiffusionTensor3D
 * \brief Represent a diffusion tensor as used in DTI images.
 *
 * This class implements a 3D symmetric tensor as it is used for representing
 * diffusion of water molecules in Diffusion Tensor Images.
 *
 * This class derive from the SymmetricSecondRankTensor and from it inherit
 * most of the Tensor related behavior. At this level we add the methods that
 * are specific to 3D and that are closely related to the concept of diffusion.
 *
 *
 * \author Jeffrey Duda from School of Engineering at University of Pennsylvania
 * \author Torsten Rohlfing from SRI International Neuroscience Program.
 *
 * This class was mostly based on files that Jeffrey Duda, Torsten Rohlfing and
 * Martin Styner contributed to the ITK users list during a discussion on
 * support for DiffusionTensorImages. A discussion on the design of this class
 * can be found in the WIKI pages of NAMIC:
 *
 * http://www.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:ITK-DiffusionTensorPixelType
 *
 * \note This work is part of the National Alliance for Medical Image
 * Computing (NAMIC), funded by the National Institutes of Health
 * through the NIH Roadmap for Medical Research, Grant U54 EB005149.
 * Information on the National Centers for Biomedical Computing
 * can be obtained from http://nihroadmap.nih.gov/bioinformatics.
 *
 *
 * \note Contributions by Torsten Rohlfing were funded by the following NIH grants
 *
 * Alcohol, HIV and the Brain,
 * NIAAA AA12999, PI: A. Pfefferbaum
 *
 * Normal Aging of Brain Structure and Function
 * NIA AG 17919, PI: E.V. Sullivan.
 *
 *
 * \par References
 * E. R. Melhem, S. Mori, G. Mukundan, M. A. Kraut, M. G. Pomper, and
 * P. C. M. van Zijl, "Diffusion tensor MR imaging of the brain and white
 * matter tractography," Am. J. Roentgenol., vol. 178, pp. 3-16, 2002.
 *
 * \sa SymmetricSecondRankTensor
 *
 * \ingroup ImageObjects   TensorObjects    Geometry
 */

template <typename TComponent>
class gtractDiffusionTensor3D : public DiffusionTensor3D<TComponent>
{
public:
  /** Standard class type alias. */
  using Self = gtractDiffusionTensor3D;
  using Superclass = itk::DiffusionTensor3D<TComponent>;
  using RealValueType = typename Superclass::RealValueType;
  using EigenValuesArrayType = typename Superclass::EigenValuesArrayType;
  /** Get the Volume Ratio from the Tensor. */
  RealValueType
  GetVolumeRatio() const;

  /** Get the value of Axial Diffusivity from the Tensor. */
  RealValueType
  GetAxialDiffusivity() const;

  /** Get the value of Radial Diffusivity from the Tensor. */
  RealValueType
  GetRadialDiffusivity() const;

  /** Get the Lattice Index from the Tensor. */
  RealValueType
  GetLatticeIndex() const;
};
} // namespace itk
#include "gtractDiffusionTensor3D.hxx"

#endif
