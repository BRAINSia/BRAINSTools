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

#ifndef __itkDtiStreamlineTrackingFilter_hxx
#define __itkDtiStreamlineTrackingFilter_hxx

#include <cmath>

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

#include "itkDtiStreamlineTrackingFilter.h"
#include "algo.h"

#include <iostream>

namespace itk
{
template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
DtiStreamlineTrackingFilter<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::DtiStreamlineTrackingFilter() : DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType,
                                                        TMaskImageType >::DtiTrackingFilterBase()
{
  this->m_CurvatureThreshold = 45;
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
void
DtiStreamlineTrackingFilter<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::Update()
{
  typedef typename Self::TensorImageType::PixelType::EigenValuesArrayType   EigenValuesArrayType;
  typedef typename Self::TensorImageType::PixelType::EigenVectorsMatrixType EigenVectorsMatrixType;

  float   anisotropy, anisotropySum(0);
  TVector vin(3), vout(3);

  typename Self::ContinuousIndexType index, tmpIndex;
  bool stop;
  bool addFiber;

  const double inRadians = this->pi / 180.0;
  double       curvatureThreshold = std::cos( this->m_CurvatureThreshold * inRadians );

  this->m_Output = vtkPolyData::New();
  this->m_TrackingDirections.clear();
  this->m_Seeds.clear();
  this->m_ScalarIP->SetInputImage(this->m_AnisotropyImage);
  this->m_VectorIP->SetInputImage(this->m_TensorImage);
  this->m_EndIP->SetInputImage(this->m_EndingRegion);
  typename Self::AnisotropyImageRegionType ImageRegion = this->m_AnisotropyImage->GetLargestPossibleRegion();

  this->m_StartIP->SetInputImage(this->m_StartingRegion);
  Self::InitializeSeeds();

  // std::cout << "Image Region for Tracking: " << ImageRegion << std::endl;
  /*** Add length and Loop Detection ***/
  while( !this->m_Seeds.empty() )
    {
    index = this->m_Seeds.back();              this->m_Seeds.pop_back();
    vin = vout  = this->m_TrackingDirections.back();  this->m_TrackingDirections.pop_back();

    stop = false;
    addFiber = false;
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

    // ////////////////////////////////////////////////////////////////////////
    // Tracking start from given 'index' and 'vout'
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

      //
      // ////////////////////////////////////////////////////////////////////////
      // evaluate the stopping criteria
      // std::cerr << "Test anisotropy " << anisotropy << " >= " <<
      // this->m_AnisotropyThreshold << std::endl;
      if( anisotropy >= this->m_AnisotropyThreshold )
        {
        if( this->m_EndIP->EvaluateAtContinuousIndex(index) >= 0.5 )
          {
          stop = true;
          addFiber = true;
          // std::cerr << "Stop - Found End" << std::endl;
          }
        // std::cerr << "Test Pathlength " << pathLength << " > " <<
        // this->m_MaximumLength << std::endl;

        if( pathLength > this->m_MaximumLength )
          {
          stop = true;
          // std::cerr << "Stop Max Length" << std::endl;
          }

        //
        // ////////////////////////////////////////////////////////////////////////
        // forward propagating
        if( currentPointId )
          {
          anisotropySum = anisotropySum + anisotropy;
          }
        else
          {
          anisotropySum = anisotropy;
          }
        fiberAnisotropy->InsertNextValue( anisotropy );
        fiberAnisotropySum->InsertNextValue( anisotropySum );

        typename Self::PointType p;
        this->ContinuousIndexToMM(index, p);
        fiber->InsertNextPoint( p.GetDataPointer() );
        currentPointId++;

        EigenValuesArrayType   eigenValues;
        EigenVectorsMatrixType eigenVectors;
        typename Self::TensorImagePixelType tensorPixel = this->m_VectorIP->EvaluateAtContinuousIndex(index);

        TMatrix fullTensorPixel(3, 3); fullTensorPixel = Tensor2Matrix( tensorPixel );
        fiberTensors->InsertNextTypedTuple( fullTensorPixel.data_block() );

        tensorPixel.ComputeEigenAnalysis(eigenValues, eigenVectors);

        //
        // ////////////////////////////////////////////////////////////////////////
        // Get major vector
        TVector e2(3); e2[0] = eigenVectors[2][0]; e2[1] = eigenVectors[2][1]; e2[2] = eigenVectors[2][2];
        if( dot_product(vin, e2) < 0 )
          {
          e2 *= -1;
          }

        //
        // ////////////////////////////////////////////////////////////////////////
        // Choose an outgoing direction
        double vin_dot_e2 = dot_product(vin, e2);
        // std::cerr << "Compare " << vin_dot_e2 << " > " <<
        // this->m_CurvatureThreshold << std::endl;
        // std::cerr << "Vin " << vin << " E2 " << e2 << std::endl;
        if( vin_dot_e2 > curvatureThreshold )
          {
          // Use TEND ???
          if( this->m_UseTend )
            {
            this->ApplyTensorDeflection(vin, fullTensorPixel, e2, vout);
            }
          else
            {
            vout = e2;
            }
          //
          // ////////////////////////////////////////////////////////////////////////
          // Calculate the new index
          this->StepIndex(tmpIndex, index, vout);
          pathLength += this->m_StepSize;

          //
          // ////////////////////////////////////////////////////////////////////////
          // Update the current index
          index = tmpIndex;
          vin = vout;
          }
        else  // Curvature Threshold
          {
          stop = true;
          // std::cerr << "Stop Curvature" << std::endl;
          }
        }
      else   // Anisotropy Threshold
        {
        stop = true;
        // std::cerr << "Stop Anisotropy" << std::endl;
        }
      }

    if( addFiber && ( pathLength >= this->m_MinimumLength ) )
      {
      std::cerr << "Fiber (" << pathLength << "); ";
      this->AddFiberToOutput( fiber, fiberTensors );
      }
    }
}
} // end namespace itk
#endif
