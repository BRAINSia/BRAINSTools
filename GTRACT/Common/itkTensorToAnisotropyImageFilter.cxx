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

#ifndef __itkTensorToAnisotropyImageFilter_hxx
#define __itkTensorToAnisotropyImageFilter_hxx

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include <itkIOCommon.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"

#include "itkTensorToAnisotropyImageFilter.h"
#include "algo.h"

#include <iostream>

namespace itk
{
TensorToAnisotropyImageFilter
::TensorToAnisotropyImageFilter()
{
  m_AnisotropyType = FRACTIONAL_ANISOTROPY;
}

void
TensorToAnisotropyImageFilter
::Update()
{
  InputImageRegionType ImageRegion = m_Input->GetLargestPossibleRegion();

  m_Output = OutputImageType::New();

  m_Output->SetRegions( ImageRegion );
  m_Output->CopyInformation( m_Input );
  m_Output->Allocate();

  switch( m_AnisotropyType )
    {
    case MEAN_DIFFUSIVITY:
      {
      computVoxelIsotropy();
      }
      break;
    case FRACTIONAL_ANISOTROPY:
    case RELATIVE_ANISOTROPY:
    case VOLUME_RATIO:
    case RADIAL_DIFFUSIVITY:
    case AXIAL_DIFFUSIVITY:
      {
      computSimpleVoxelAnisotropy();
      }
      break;
    case COHERENCE_INDEX:
    case LATTICE_INDEX:
      {
      computNeighborhoodVoxelAnisotropy();
      }
      break;
    default:
      {
      }
      break;
    }

  // Set Meta Data Orientation Information
  m_Output->SetMetaDataDictionary( m_Input->GetMetaDataDictionary() );
}

void
TensorToAnisotropyImageFilter
::computVoxelIsotropy()
{
  typedef itk::ImageRegionIteratorWithIndex<OutputImageType> IteratorType;
  IteratorType it( m_Output, m_Output->GetLargestPossibleRegion() );

  OutputImageIndexType index;
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    InputPixelType currentVoxel = m_Input->GetPixel( index );
    float          adc = 0;

    if( currentVoxel.GetNorm() != 0 )
      {
      adc = ( currentVoxel[0] + currentVoxel[1] + currentVoxel[2] ) / 3.0;
      }
    it.Set( adc );
    }
}

void
TensorToAnisotropyImageFilter
::computSimpleVoxelAnisotropy()
{
  typedef itk::ImageRegionIteratorWithIndex<OutputImageType> IteratorType;
  IteratorType it( m_Output, m_Output->GetLargestPossibleRegion() );

  OutputImageIndexType index;
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    InputPixelType currentVoxel = m_Input->GetPixel( index );
    float          fa = 0;

    if( currentVoxel.GetNorm() != 0 )
      {
      TVector eig = Eigen_Value( Tensor2Matrix( currentVoxel ) );

      switch( m_AnisotropyType )
        {
        case FRACTIONAL_ANISOTROPY:
          {
          fa = FA(eig);
          }
          break;
        case RELATIVE_ANISOTROPY:
          {
          fa = RA(eig);
          }
          break;
        case VOLUME_RATIO:
          {
          fa = VR(eig);
          }
          break;
        case AXIAL_DIFFUSIVITY:
          {
          fa = AxialDiffusivity(eig);
          }
          break;
        case RADIAL_DIFFUSIVITY:
          {
          fa = RadialDiffusivity(eig);
          }
          break;
        default:
          {
          fa = 0;
          }
          break;
        }
      }
    it.Set(fa);
    }
}

void
TensorToAnisotropyImageFilter
::computNeighborhoodVoxelAnisotropy()
{
  m_Output->FillBuffer(0.0);
  typedef itk::ConstNeighborhoodIterator<InputImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType it;

  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);
  radius[2] = 0;

  // boundary condition
  typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> FaceCalculatorType;
  FaceCalculatorType               faceCalculator;
  FaceCalculatorType::FaceListType faceList;
  faceList = faceCalculator(m_Input, m_Input->GetLargestPossibleRegion(), radius);
  FaceCalculatorType::FaceListType::iterator fit;
  for( fit = faceList.begin(); fit != faceList.end(); ++fit )
    { // This is temporary, further consideration on boundary condition needed
    it = NeighborhoodIteratorType( radius, m_Input, *fit );
    for( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      TVector center = it.GetCenterPixel().GetVnlVector();
      float   ai = 0;
      // ////////////////////////////////////////////////////////////////////////
      if( !center.is_zero() )
        {
        float   sum = 0;
        float   coef = 0;
        TVector neighbor;
        for( int i = 0; i <= 8; i++ )
          {
          if( i == 4 )
            {
            continue;
            }

          neighbor = it.GetPixel(i).GetVnlVector();
          if( !neighbor.is_zero() )
            {
            float temp;
            float a = 1;
            if( ( i % 2 ) == 0 )
              {
              a = 0.7071;
              }

            switch( m_AnisotropyType )
              {
              case COHERENCE_INDEX:
                {
                temp = CI(center, neighbor);
                }
                break;
              case LATTICE_INDEX:
                {
                temp = LI(center, neighbor);
                }
                break;
              default:
                {
                temp = 0;
                }
                break;
              }

            sum += a * temp;
            coef += a;
            }
          }

        // Cut off the value that < 0, It's right or wrong?
        if( ( coef != 0 ) & ( sum > 0 ) )
          {
          ai = sum / coef;
          }
        }
      // ////////////////////////////////////////////////////////////////////////
      m_Output->SetPixel(it.GetIndex(), ai);
      }
    }
}
} // end namespace itk
#endif
