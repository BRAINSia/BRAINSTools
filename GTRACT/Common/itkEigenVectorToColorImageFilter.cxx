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

#ifndef __itkEigenVectorToColorImageFilter_hxx
#define __itkEigenVectorToColorImageFilter_hxx

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include <itkIOCommon.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"

#include "itkEigenVectorToColorImageFilter.h"
#include "algo.h"

#include <iostream>

namespace itk
{
EigenVectorToColorImageFilter
::EigenVectorToColorImageFilter()
{
  m_TensorShapeType = PRIMARY_EIGENVECTOR;
}

void
EigenVectorToColorImageFilter
::Update()
{
  InputImageType::RegionType ImageRegion = m_Input->GetLargestPossibleRegion();

  m_Output = OutputImageType::New();
  m_Output->SetRegions( ImageRegion );
  m_Output->CopyInformation( m_Input );
  m_Output->Allocate();

  typedef itk::ImageRegionIteratorWithIndex<OutputImageType> IteratorType;

  IteratorType it( m_Output, ImageRegion );

  OutputImageType::IndexType index;
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    InputPixelType  tensor = m_Input->GetPixel(index);
    OutputPixelType currentVoxel;
    currentVoxel.Fill(0);
    if( tensor.GetNorm() != 0 )
      {
      TMatrix                          M = Tensor2Matrix(tensor);
      vnl_symmetric_eigensystem<float> eig( M);
      TVector                          e(3);
      float                            ai = FA( Eigen_Value(M) );

      switch( m_TensorShapeType )
        {
        case PRIMARY_EIGENVECTOR:
        case SECONDARY_EIGENVECTOR:
        case TERTIARY_EIGENVECTOR:
          {
          e = eig.get_eigenvector( static_cast<int>( m_TensorShapeType ) );
          }
          break;
        case TENSOR_SHAPE:
          TVector eigV(3);
          eigV(2) = eig.get_eigenvalue(2);
          eigV(1) = eig.get_eigenvalue(1);
          eigV(0) = eig.get_eigenvalue(0);
          e = TensorShape(eigV);
        }
      e *= ai;

      /**************************************************************
      The standard orientation internal is Coronal: The orientation is thus
        Red: Right - Left - X axis in Image
        Green: Anterior-Posterior - z axis in Image
        Blue: Superior - Inferior - y axis in Image
      **************************************************************/
      currentVoxel.SetRed( static_cast<unsigned char>( 255 * vcl_sqrt( vcl_fabs( e(0) ) ) ) );
      currentVoxel.SetGreen( static_cast<unsigned char>( 255 * vcl_sqrt( vcl_fabs( e(2) ) ) ) );
      currentVoxel.SetBlue( static_cast<unsigned char>( 255 * vcl_sqrt( vcl_fabs( e(1) ) ) ) );
      currentVoxel.SetAlpha(255);
      }
    it.Set(currentVoxel);
    }

  // Set Meta Data Orientation Information
  m_Output->SetMetaDataDictionary( m_Input->GetMetaDataDictionary() );
}
} // end namespace itk
#endif
