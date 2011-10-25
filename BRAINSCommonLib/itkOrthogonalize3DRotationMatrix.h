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
