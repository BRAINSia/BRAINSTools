#include "itkOrthogonalize3DRotationMatrix.h"

namespace itk
{
itk::Matrix<double, 3, 3> Orthogonalize3DRotationMatrix( const itk::Matrix<double, 3, 3> & rotator)
{
  vnl_svd<double>                             decomposition( rotator.GetVnlMatrix(), -1e-6 );
  vnl_diag_matrix<vnl_svd<double>::singval_t> Winverse( decomposition.Winverse() );

  vnl_matrix<double> W(3, 3);
  W.fill( double(0) );
  for( unsigned int i = 0; i < 3; ++i )
    {
    if( decomposition.Winverse() (i, i) != 0.0 )
      {
      W(i, i) = 1.0;
      }
    }

  const vnl_matrix<double> result(
    decomposition.U() * W * decomposition.V().conjugate_transpose() );

  //    std::cout << " svd Orthonormalized Rotation: " << std::endl
  //      << result << std::endl;
  itk::Matrix<double, 3, 3> Orthog;
  Orthog.operator=(result);

  return Orthog;
}
}
