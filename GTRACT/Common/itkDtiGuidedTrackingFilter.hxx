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

#ifndef __itkDtiGuidedTrackingFilter_hxx
#define __itkDtiGuidedTrackingFilter_hxx

#include "vtkFloatArray.h"
#include "vtkPolyLine.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include <itkIOCommon.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"

#include "itkDtiGuidedTrackingFilter.h"
#include "algo.h"

#include <iostream>

namespace itk
{
template <typename TTensorImageType, typename TAnisotropyImageType, typename TMaskImageType>
DtiGuidedTrackingFilter<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::DtiGuidedTrackingFilter() : DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType,
                                                    TMaskImageType >::DtiTrackingFilterBase()
{
  this->m_CurvatureThreshold       = 45.0;
  this->m_GuidedCurvatureThreshold = 30.0;
  this->m_MaximumGuideDistance     = 12.0;
}

template <typename TTensorImageType, typename TAnisotropyImageType, typename TMaskImageType>
void
DtiGuidedTrackingFilter<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::SetGuideFiber( GuideFiberType guideFiber)
{
  this->m_GuideFiber = guideFiber;
}

template <typename TTensorImageType, typename TAnisotropyImageType, typename TMaskImageType>
void
DtiGuidedTrackingFilter<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::Update()
{
  using EigenValuesArrayType = typename Self::TensorImageType::PixelType::EigenValuesArrayType;
  using EigenVectorsMatrixType = typename Self::TensorImageType::PixelType::EigenVectorsMatrixType;

  // std::cout << this->m_AnisotropyImage;

  this->m_Output = vtkPolyData::New();
  this->m_ScalarIP->SetInputImage(this->m_AnisotropyImage);
  this->m_VectorIP->SetInputImage(this->m_TensorImage);
  this->m_EndIP->SetInputImage(this->m_EndingRegion);
  typename Self::AnisotropyImageRegionType ImageRegion = this->m_AnisotropyImage->GetLargestPossibleRegion();

  this->m_StartIP->SetInputImage(this->m_StartingRegion);
  Self::InitializeSeeds();

  // ////////////////////////////////////////////////////////////////////////
  // Initialize some parameters
  float   anisotropy, anisotropySum(0);
  TVector vin(3), vout(3), vguide(3);

  typename Self::ContinuousIndexType index, tmpIndex;
  bool stop;
  int  count = 0;

  const double inRadians = this->pi / 180.0;
  double       curvatureThreshold = std::cos( this->m_CurvatureThreshold * inRadians );
  double       guidedCurvatureThreshold = std::cos( this->m_GuidedCurvatureThreshold * inRadians );

  // std::cout << "#Seeds = " << this->m_Seeds.size() << std::endl;

  // ////////////////////////////////////////////////////////////////////////
  // For each seed point, start guided tracking
  while( !this->m_Seeds.empty() )
    {
    index = this->m_Seeds.back();              this->m_Seeds.pop_back();
    vin = vout  = this->m_TrackingDirections.back();  this->m_TrackingDirections.pop_back();

    stop = false;
    // addFiber = false;
    vtkPoints *    fiber = vtkPoints::New();
    vtkFloatArray *fiberTensors = vtkFloatArray::New();
    fiberTensors->SetName("Tensors");
    fiberTensors->SetNumberOfComponents(9);
    vtkFloatArray *fiberAnisotropy = vtkFloatArray::New();
    fiberAnisotropy->SetName("Anisotropy");
    vtkFloatArray *fiberAnisotropySum = vtkFloatArray::New();
    fiberAnisotropySum->SetName("Anisotropy-Sum");
    int   currentPointId = 0;
    float pathLength = 0.0;

    typename Self::PointType p2;
    double p1[3];
    this->m_GuideFiber->GetPoint(0, p1);
    this->ContinuousIndexToMM( index, p2 );
    typename Self::ContinuousIndexType index1;
    this->MMToContinuousIndex( p1, index1 );

    // std::cout << "Guide Index: " << index1 << std::endl;
    // std::cout << "Guide Point 0: " << p1[0] << " " << p1[1] << " " << p1[2]
    // << std::endl;
    // std::cout << "Seed Index: " << index << std::endl;
    // std::cout << "Seed Point: " << p2 << std::endl;

    /***VAM - MaxDistance is now defined by the user */
    // float MaxDist =
    //
    // std::sqrt(pow((double)(p1[0]-p2[0]),2.0)+pow((double)(p1[1]-p2[1]),2.0)+pow((double)(p1[2]-p2[2]),2.0));
    // MaxDist *= 1.5;
    double MaxDist = this->m_MaximumGuideDistance;
    // std::cout << "Max Distance: " << MaxDist << std::endl;

    while( !stop )
      {
      if( ImageRegion.IsInside(index) )
        {
        anisotropy = this->m_ScalarIP->EvaluateAtContinuousIndex(index);
        }
      else
        {
        anisotropy = -1;
        }

      //
      // ////////////////////////////////////////////////////////////////////////
      // Evaluate the stopping criteria: is below fa threshold? is outside image
      // region?
      if( anisotropy >= this->m_AnisotropyThreshold )
        {
        if( currentPointId == 0 )
          {
          anisotropySum = anisotropy;
          }
        else
          {
          anisotropySum += anisotropy;
          }

        typename Self::PointType p;
        this->ContinuousIndexToMM( index, p );
        fiber->InsertNextPoint( p.GetDataPointer() );

        currentPointId++;
        fiberAnisotropy->InsertNextValue( anisotropy );
        fiberAnisotropySum->InsertNextValue( anisotropySum );

        // std::cout << "\tFiber Point: " << index << std::endl;

        //
        // ////////////////////////////////////////////////////////////////////////
        // Seeking guidance
        bool isGuided = GuideDirection(index, this->m_GuideFiber, MaxDist, vguide);

        EigenValuesArrayType   eigenValues;
        EigenVectorsMatrixType eigenVectors;
        typename Self::TensorImagePixelType tensorPixel = this->m_VectorIP->EvaluateAtContinuousIndex(index);

        TMatrix fullTensorPixel(3, 3); fullTensorPixel = Tensor2Matrix( tensorPixel );
        fiberTensors->InsertNextTypedTuple( fullTensorPixel.data_block() );

        tensorPixel.ComputeEigenAnalysis(eigenValues, eigenVectors);

        TVector e2(3); e2[0] = eigenVectors[2][0]; e2[1] = eigenVectors[2][1]; e2[2] = eigenVectors[2][2];
        // std::cout << "\tEigen Vector " << e2 << " Guide Direction " << vguide
        // << std::endl;
        if( isGuided )
          {
          // std::cout << "\tGuided Fiber: " << std::endl;

          if( dot_product(e2, vin) < 0 )
            {
            e2 *= -1;
            }

          if( dot_product(vguide, vin) < 0 )
            {
            vguide *= -1;
            }

          if( dot_product(e2, vguide) < guidedCurvatureThreshold )
            {
            vout = vguide; // using guiding direction
            // std::cout << "\tUsing Guide Direction: " << vguide << std::endl;
            }
          else
            {
            // std::cout << "\tUsing EigenVector Direction: " << e2 <<
            // std::endl;
            // Use tend???
            if( this->m_UseTend )
              {
              this->ApplyTensorDeflection(vin, fullTensorPixel, e2, vout);
              }
            else
              {
              vout  = e2;
              }
            }

          // std::cout << "\tOut Direction: " << vout << std::endl;
          //
          // ////////////////////////////////////////////////////////////////////////
          // Update Index
          this->StepIndex(tmpIndex, index, vout);
          pathLength += this->m_StepSize;
          index = tmpIndex;
          vin = vout;
          // std::cout << "New Index: " << index << std::endl;
          }
        else
          {
          //
          // ////////////////////////////////////////////////////////////////////////
          // Unguided -- can't use the guide
          // std::cout << "\tUnguided Fiber: " << std::endl;
          //
          // ////////////////////////////////////////////////////////////////////////
          // Get the principle eigen vector at the current point

          if( dot_product(vin, e2) < 0 )
            {
            e2 *= -1;
            }

          // Check the Curvature Threshold
          if( dot_product(vin, e2) < curvatureThreshold )
            {
            if( this->m_UseTend )
              {
              this->ApplyTensorDeflection(vin, fullTensorPixel, e2, vout);
              }
            else
              {
              vout = e2;
              }

            this->StepIndex(tmpIndex, index, vout);
            pathLength += this->m_StepSize;
            index = tmpIndex;
            vin = vout;
            }
          else
            {
            // std::cout << "Abandon Unguided Fiber Below Curvature Threshold"
            // << std::endl;
            stop = true;
            }
          }
        //
        // ////////////////////////////////////////////////////////////////////////
        }
      else
        {
        stop = true;
        // std::cout << "Abandon Fiber Below Anisotropy Threshold" << std::endl;
        }

      if( ( this->m_EndIP->EvaluateAtContinuousIndex(index) >= 0.5 ) && ( pathLength >= this->m_MinimumLength ) )
        {
        std::cerr << "\tFiber-" << count << " (" << currentPointId <<  "); ";
        this->AddFiberToOutput( fiber, fiberTensors );
        stop = true;
        count++;
        }

      // Check for loops if selected by the user
      if( this->m_UseLoopDetection )
        {
        if( Self::IsLoop(fiber) )
          {
          stop = true;
          }
        }

      // Check fiber length
      if( pathLength > this->m_MaximumLength )
        {
        // std::cout << "Abandon Max Length" << fiber->GetNumberOfPoints() <<
        // std::endl;
        stop = true;
        }
      } // Fiber Path Loop
    }   // Seed Loop
}

template <typename TTensorImageType, typename TAnisotropyImageType, typename TMaskImageType>
bool
DtiGuidedTrackingFilter<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::GuideDirection(typename Self::ContinuousIndexType index,
                 GuideFiberType centerFiber,
                 const float MaxDist,
                 TVector & vguide)
{
  TVector direction(3);

  float minDist = MaxDist;

  typename Self::PointType p1, p2, p3;
  // this->m_AnisotropyImage->TransformContinuousIndexToPhysicalPoint(index,p2);
  // std::cout << "Current Point " << index << std::endl;
  for( int i = 0; i < centerFiber->GetNumberOfPoints(); i++ )
    {
    double currentPoint[3];
    centerFiber->GetPoint(i, currentPoint);
    p1[0] = currentPoint[0]; p1[1] = currentPoint[1]; p1[2] = currentPoint[2];
    typename Self::ContinuousIndexType index1;
    this->m_AnisotropyImage->TransformPhysicalPointToContinuousIndex(p1, index1);

    float dist = std::sqrt( std::pow( (double)( index1[0] - index[0] ), 2.0 )
                           + std::pow( (double)( index1[1] - index[1] ), 2.0 )
                           + std::pow( (double)( index1[2] - index[2] ), 2.0 ) );
    if( dist < minDist )
      {
      minDist = dist;
      typename Self::ContinuousIndexType index3;
      if( i == centerFiber->GetNumberOfPoints() - 1 )
        {
        centerFiber->GetPoint(i - 1, currentPoint);
        p3[0] = currentPoint[0]; p3[1] = currentPoint[1]; p3[2] = currentPoint[2];
        this->m_AnisotropyImage->TransformPhysicalPointToContinuousIndex(p3, index3);
        for( int j = 0; j < 3; j++ )
          {
          direction[j] = index1[j] - index3[j];
          }
        }
      else
        {
        centerFiber->GetPoint(i + 1, currentPoint);
        p3[0] = currentPoint[0]; p3[1] = currentPoint[1]; p3[2] = currentPoint[2];
        this->m_AnisotropyImage->TransformPhysicalPointToContinuousIndex(p3, index3);
        for( int j = 0; j < 3; j++ )
          {
          direction[j] = index3[j] - index1[j];
          }
        }
      }
    }
  // std::cout << "Current Point " << index << " Min dist : " << minDist << "
  // Max Dist: " << MaxDist << std::endl;

  if( minDist >= MaxDist )
    {
    return false;
    }
  else
    {
    direction.normalize();
    vguide = direction;
    return true;
    }
}
} // end namespace itk
#endif
