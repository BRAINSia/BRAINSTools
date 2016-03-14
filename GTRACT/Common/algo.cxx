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

#include "algo.h"

TMatrix Matrix_Inverse( TMatrix M )
{
  //  const int NumberOfDirections=M.rows();
  TMatrix M_T = M.transpose();

  vnl_matrix_inverse<float> M_Inverse(M_T * M);
  return M_Inverse.inverse() * M_T;
}

float My_lsf( TVector x, TVector y)
{
  unsigned int N = x.size();

  if( N != y.size() )
    {
    throw;
    }
  float a, b, c, d;
  a = b = c = d = 0;
  for( unsigned int i = 0; i < N; i++ )
    {
    a += x(i) * y(i);
    b += x(i) * x(i);
    c += x(i);
    d += y(i);
    }
  return ( a / c - d / N ) / ( b / c - c / N );
}

TVector Eigen_Value( TMatrix M )
{
  vnl_symmetric_eigensystem<float> eig(M);
  TVector                          result(3);
  for( unsigned int i = 0; i < 3; i++ )
    {
    result(i) = eig.get_eigenvalue(i);
    }
  return result;
}

TMatrix Tensor2Matrix(TVector ADCe)
{
  TMatrix M(3, 3);

  M(0, 0) = ADCe(0);  M(0, 1) = ADCe(1);  M(0, 2) = ADCe(2);
  M(1, 0) = ADCe(1);  M(1, 1) = ADCe(3);  M(1, 2) = ADCe(4);
  M(2, 0) = ADCe(2);  M(2, 1) = ADCe(4);  M(2, 2) = ADCe(5);
  return M;
}

TMatrix Tensor2Matrix(itk::FixedArray<float, 6> ADCe)
{
  TMatrix M(3, 3);

  M(0, 0) = ADCe[0];  M(0, 1) = ADCe[1];  M(0, 2) = ADCe[2];
  M(1, 0) = ADCe[1];  M(1, 1) = ADCe[3];  M(1, 2) = ADCe[4];
  M(2, 0) = ADCe[2];  M(2, 1) = ADCe[4];  M(2, 2) = ADCe[5];
  return M;
}

TMatrix Tensor2Matrix(itk::FixedArray<double, 6> ADCe)
{
  TMatrix M(3, 3);

  M(0, 0) = ADCe[0];  M(0, 1) = ADCe[1];  M(0, 2) = ADCe[2];
  M(1, 0) = ADCe[1];  M(1, 1) = ADCe[3];  M(1, 2) = ADCe[4];
  M(2, 0) = ADCe[2];  M(2, 1) = ADCe[4];  M(2, 2) = ADCe[5];
  return M;
}

TVector DD(TVector ADC1, TVector ADC2)
{
  TVector temp(2, 0);

  vnl_symmetric_eigensystem<float> eig1( Tensor2Matrix(ADC1) );
  vnl_symmetric_eigensystem<float> eig2( Tensor2Matrix(ADC2) );

  TVector V1(3);
  for( int i = 0; i < 3; i++ )
    {
    V1(i) = eig1.get_eigenvalue(i);
    }
  TVector V2(3);
  for( int i = 0; i < 3; i++ )
    {
    V2(i) = eig2.get_eigenvalue(i);
    }
  for( int i = 0; i < 3; i++ )
    {
    for( int j = 0; j < 3; j++ )
      {
      temp(0) += V1(i) * V2(j)
        * std::pow(dot_product( eig1.get_eigenvector(i), eig2.get_eigenvector(j) ), 2);
      }
    }

  temp(1) = V1.sum() * V2.sum() / 3;

  return temp;
}

float CI(TVector ADC1, TVector ADC2)
{
  TVector     temp = DD(ADC1, ADC2);
  const float result = ( temp(0) - temp(1) ) / temp(0);

  return result;
}

float LI(TVector ADC1, TVector ADC2)
{
  TVector t0 = DD(ADC1, ADC2);

  // When two tensor is equal, then the DD operation could be simplified
  vnl_symmetric_eigensystem<float> eig1( Tensor2Matrix(ADC1) );
  const float                      t1 =
    std::pow(eig1.get_eigenvalue(0), 2) + std::pow(eig1.get_eigenvalue(1), 2) + std::pow( eig1.get_eigenvalue(2), 2);
  const vnl_symmetric_eigensystem<float> eig2( Tensor2Matrix(ADC2) );
  float                                  t2 =
    std::pow(eig2.get_eigenvalue(0), 2) + std::pow(eig2.get_eigenvalue(1), 2) + std::pow( eig2.get_eigenvalue(2), 2);

  float result;
  if( ( ( t0(0) - t0(1) ) < 0 ) || ( t0(0) < 0 ) )
    {
    result = 0;
    }
  else
    {
    result = 0.612372 * std::sqrt( t0(0) - t0(1) ) / std::sqrt( t0(0) )
      + 0.75 * ( t0(0) - t0(1) ) / ( std::sqrt(t1 * t2) );
    }
  return result;
}

float FA( TVector eig )
{
  const float mean = ( eig(0) + eig(1) + eig(2) ) / 3;
  const float up = std::sqrt( std::pow(eig(0) - mean, 2) + std::pow(eig(1) - mean, 2) + std::pow(eig(2) - mean, 2) );
  // float up = std::sqrt( std::pow(eig(0)-eig(1), 2) + std::pow(eig(1)-eig(2),2) +
  // std::pow(eig(2)-eig(0),2) );
  const float bottom = std::sqrt( std::pow(eig(0), 2) + std::pow(eig(1), 2) + std::pow(eig(2), 2) );
  float       result = std::sqrt(1.5) * up / bottom;

  if( result > 1.0 )
    {
    result = 1.0;
    }
  return result;
}

float RA( TVector eig )
{
  const float mean = ( eig(0) + eig(1) + eig(2) ) / 3;
  const float up = std::sqrt( std::pow(eig(0) - mean, 2) + std::pow(eig(1) - mean, 2) + std::pow(eig(2) - mean, 2) );
  const float bottom = std::sqrt( 3.0) * mean;

  return up / bottom;
}

float VR( TVector eig )
{
  const float mean = ( eig(0) + eig(1) + eig(2) ) / 3;
  const float up = eig(0) * eig(1) * eig(2);
  const float bottom = std::pow( mean, 3.0F );

  return up / bottom;
}

TVector TensorShape(TVector eigV)
{
  TVector result(3);

  result(0) = ( eigV(2) - eigV(1) ) / eigV.sum();
  result(1) = 2 * ( eigV(1) - eigV(0) ) / eigV.sum();
  result(2) = 3 * eigV(0) / eigV.sum();
  return result;
}

float MeanDiffusivity( TVector eig )
{
  const float mean = ( eig(0) + eig(1) + eig(2) ) / 3;

  return mean;
}

float AxialDiffusivity( TVector eig )
{
  const float result = eig(2);

  return result;
}

float RadialDiffusivity( TVector eig )
{
  const float result = ( eig(0) + eig(1) ) / 2.0;

  return result;
}
