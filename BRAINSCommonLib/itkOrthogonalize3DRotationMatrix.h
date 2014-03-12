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
#ifndef __itkOrthogonalize3DRotationMatrix_h__
#define __itkOrthogonalize3DRotationMatrix_h__
#include "itkMatrix.h"

namespace itk
{
/**
 * Orthogonalize3DRotationMatrix
 *
 * This function will take in a rotation matrix
 * and generate a version that is orthogonal
 * to a very high precision while minimizing the
 * actual change in rotational values.
 */
itk::Matrix<double, 3, 3> Orthogonalize3DRotationMatrix( const itk::Matrix<double, 3, 3> & rotator);
}

#endif // __itkOrthogonalize3DRotationMatrix_h__
