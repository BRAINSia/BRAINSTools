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
#include "itkOrthogonalize3DRotationMatrix.h"

namespace itk
{
itk::Matrix<double, 3, 3>
Orthogonalize3DRotationMatrix(const itk::Matrix<double, 3, 3> & rotator)
{
  vnl_svd<double>                                   decomposition(rotator.GetVnlMatrix().as_matrix(), -1e-6);
  const vnl_diag_matrix<vnl_svd<double>::singval_t> Winverse(decomposition.Winverse());

  vnl_matrix<double> W(3, 3);
  W.fill(static_cast<double>(0));
  for (unsigned int i = 0; i < 3; ++i)
  {
    if (decomposition.Winverse()(i, i) != 0.0)
    {
      W(i, i) = 1.0;
    }
  }

  const vnl_matrix<double> result(decomposition.U() * W * decomposition.V().conjugate_transpose());

  //    std::cout << " svd Orthonormalized Rotation: " << std::endl
  //      << result << std::endl;
  itk::Matrix<double, 3, 3> Orthog;
  Orthog.operator=(result);

  return Orthog;
}
} // namespace itk
