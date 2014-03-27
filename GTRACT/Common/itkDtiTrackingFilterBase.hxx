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

#ifndef __itkDtiTrackingFilterBase_hxx
#define __itkDtiTrackingFilterBase_hxx

#include "vtkAppendPolyData.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkVersionMacros.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
// #include <itkIOCommon.h>
// #include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"

#include "itkDtiTrackingFilterBase.h"
// #include "algo.h"


#include <iostream>

namespace itk
{
template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::DtiTrackingFilterBase()
{
  m_UseTend = false;
  m_UseLoopDetection = true;
  m_TendG = 1.0;
  m_TendF = 0.0;
  m_StepSize = 1.0;
  m_MaximumLength = 100.0;
  m_MinimumLength = 0.0;
  m_AnisotropyThreshold = 0.3;
  m_SeedThreshold = 0.5;
  m_ScalarIP    = ScalarIPType::New();
  m_VectorIP    = VectorIPType::New();
  m_StartIP             = Self::MaskIPType::New();
  m_EndIP               = Self::MaskIPType::New();
  pi = 3.14159265358979323846;
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
void
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::ContinuousIndexToMM(typename Self::ContinuousIndexType & index, PointType & p)
{
  this->m_AnisotropyImage->TransformContinuousIndexToPhysicalPoint(index, p);
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
void
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::MMToContinuousIndex(PointType & p, typename Self::ContinuousIndexType & index)
{
  this->m_AnisotropyImage->TransformPhysicalPointToContinuousIndex(p, index);
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
void
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::MMToContinuousIndex(double *pt, typename Self::ContinuousIndexType & index)
{
  PointType p; p[0] = pt[0]; p[1] = pt[1]; p[2] = pt[2];

  this->m_AnisotropyImage->TransformPhysicalPointToContinuousIndex(p, index);
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
void
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::StepIndexInPointSpace(typename Self::ContinuousIndexType & newIndex,
                        typename Self::ContinuousIndexType & oldIndex,
                        TVector & vec)
{
  PointType oldpt, newpt;

  this->m_AnisotropyImage->TransformContinuousIndexToPhysicalPoint(oldIndex, oldpt);
  // std::cerr << "Converted " << oldIndex << " to " << oldpt << std::endl;
  // Calculate the new point
  for( int i = 0; i < 3; i++ )
    {
    newpt[i] = oldpt[i] + vec[i] * this->m_StepSize;
    }

  this->m_AnisotropyImage->TransformPhysicalPointToContinuousIndex(newpt, newIndex);
  // std::cerr << "Converted " << newpt << " to " << newIndex << std::endl;
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
void
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::StepIndex(typename Self::ContinuousIndexType & newIndex,
            typename Self::ContinuousIndexType & oldIndex,
            TVector & vec)
{
  typename Self::AnisotropyImageType::SpacingType spacing = this->m_AnisotropyImage->GetSpacing();
  // Calculate the new index
  for( int i = 0; i < 3; i++ )
    {
    newIndex[i] = oldIndex[i] + vec[i] * this->m_StepSize / spacing[i];
    }
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
void
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::ApplyTensorDeflection(TVector & vin, TMatrix & fullTensorPixel, TVector & e2, TVector & vout)
{
  TVector deflection(3); deflection = fullTensorPixel * vin;

  deflection.normalize();

  vout = e2 * this->m_TendF + ( vin * ( 1 - this->m_TendG ) + deflection * this->m_TendG ) * ( 1 - this->m_TendF );
  vout.normalize();
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
typename DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>::DtiFiberType
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::GetOutput()
{
  return m_Output;
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
bool
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::IsLoop(vtkPoints *fiber, double tolerance)
{
  double p1[3], p2[3];

  const double tol2 = tolerance * tolerance;
  const int    numPts = fiber->GetNumberOfPoints();

  fiber->GetPoint(numPts - 1, p1);
  for( int i = numPts - 2; i >= 0; i-- )
    {
    fiber->GetPoint(i, p2);
    const double distance
      = ( p1[0]
          - p2[0] ) * ( p1[0] - p2[0] ) + ( p1[1] - p2[1] ) * ( p1[1] - p2[1] ) + ( p1[2] - p2[2] ) * ( p1[2] - p2[2] );
    if( distance < tol2 )
      {
      return true;
      }
    }
  return false;
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
void
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::InitializeSeeds()
{
  // ////////////////////////////////////////////////////////////////////////
  // Initialize the seed points

  typedef typename Self::TensorImageType::PixelType::EigenValuesArrayType   EigenValuesArrayType;
  typedef typename Self::TensorImageType::PixelType::EigenVectorsMatrixType EigenVectorsMatrixType;

  typedef itk::ImageRegionConstIterator<MaskImageType> ConstMaskIteratorType;
  ConstMaskIteratorType maskIt( m_StartingRegion, m_StartingRegion->GetLargestPossibleRegion() );

  int count = 0;
  int maskcount = 0;
  for( maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt )
    {
    typename Self::ContinuousIndexType        seed;
    typename ConstMaskIteratorType::IndexType pos = maskIt.GetIndex();
    seed[0] = pos[0];  seed[1] = pos[1];  seed[2] = pos[2];
    const float ai = m_ScalarIP->EvaluateAtContinuousIndex(seed);
    // const float roi = m_StartIP->EvaluateAtContinuousIndex(seed);
    if( maskIt.Get() )
      {
      maskcount++;
      if( ai >= m_SeedThreshold )
        {
        EigenValuesArrayType   eigenValues;
        EigenVectorsMatrixType eigenVectors;
        typename Self::TensorImagePixelType tensorPixel = this->m_VectorIP->EvaluateAtContinuousIndex(seed);
        tensorPixel.ComputeEigenAnalysis(eigenValues, eigenVectors);
        TVector direction(3); direction[0] = eigenVectors[2][0]; direction[1] = eigenVectors[2][1]; direction[2]
          = eigenVectors[2][2];
        m_Seeds.push_back(seed);
        m_TrackingDirections.push_back(direction);
        direction *= -1;
        m_Seeds.push_back(seed);
        m_TrackingDirections.push_back(direction);
        count++;
        }
      }
    }
  std::cerr << "Number of voxels in mask: " << maskcount << std::endl;
  std::cerr << "Number of Seeds: " << count << std::endl;
}

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
void
DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
::AddFiberToOutput( vtkPoints *currentFiber, vtkFloatArray *fiberTensors )
{
  // std::cerr << "NumPts " << currentFiber->GetNumberOfPoints() << ".  ";

  vtkCellArray *line = vtkCellArray::New();

  line->InsertNextCell( currentFiber->GetNumberOfPoints() );
  for( int i = 0; i < currentFiber->GetNumberOfPoints(); i++ )
    {
    line->InsertCellPoint(i);
    }
  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints( currentFiber );
  data->SetLines(line);
  // data->GetPointData()->SetScalars(fiberAnisotropy);
  data->GetPointData()->SetTensors(fiberTensors);

  vtkAppendPolyData *append = vtkAppendPolyData::New();
#if (VTK_MAJOR_VERSION < 6)
  append->AddInput( this->m_Output );
  append->AddInput( data );
#else
  append->AddInputData( this->m_Output );
  append->AddInputData( data );
#endif
  append->Update();
  // need to erase the old m_Output
  //  vtkPolyData *former = this->m_Output;
  this->m_Output = append->GetOutput();
  // Explicit delete because they are not ITK SmartPointers.
  //  former->Delete();
  //  append->Delete();
  //  data->Delete();
  //  line->Delete();
}
} // end namespace itk
#endif
