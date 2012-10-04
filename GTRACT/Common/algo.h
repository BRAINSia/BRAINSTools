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

/** \
 *   This file implement some algorithms
 *  Like Matrix Inverse, Least Square Fit, the Calculation of the Eigen Value and FA value
 */
#ifndef __algo_H_
#define __algo_H_

#include <cmath>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_transpose.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <itkFixedArray.h>

#include "GtractTypes.h"
#include "gtractCommonWin32.h"

extern GTRACT_COMMON_EXPORT TMatrix Matrix_Inverse( TMatrix M );

extern GTRACT_COMMON_EXPORT float My_lsf( TVector x, TVector y);

extern GTRACT_COMMON_EXPORT TVector Eigen_Value( TMatrix M );

extern GTRACT_COMMON_EXPORT TMatrix Tensor2Matrix(TVector ADCe);

extern GTRACT_COMMON_EXPORT TMatrix Tensor2Matrix(itk::FixedArray<float, 6> ADCe);

extern GTRACT_COMMON_EXPORT TMatrix Tensor2Matrix(itk::FixedArray<double, 6> ADCe);

extern GTRACT_COMMON_EXPORT float FA( TVector eig );

extern GTRACT_COMMON_EXPORT float RA( TVector eig );

extern GTRACT_COMMON_EXPORT float VR( TVector eig );

extern GTRACT_COMMON_EXPORT TVector DD(TVector ADC1, TVector ADC2);

extern GTRACT_COMMON_EXPORT float CI(TVector ADC1, TVector ADC2);

extern GTRACT_COMMON_EXPORT float LI(TVector ADC1, TVector ADC2);

extern GTRACT_COMMON_EXPORT TVector TensorShape(TVector eigV);

extern GTRACT_COMMON_EXPORT float MeanDiffusivity( TVector eig );

extern GTRACT_COMMON_EXPORT float AxialDiffusivity( TVector eig );

extern GTRACT_COMMON_EXPORT float RadialDiffusivity( TVector eig );

#endif
