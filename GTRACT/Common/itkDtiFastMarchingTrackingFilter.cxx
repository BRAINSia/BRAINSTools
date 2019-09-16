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

#ifndef __itkDtiFastMarchingTrackingFilter_cxx
#define __itkDtiFastMarchingTrackingFilter_cxx

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodAlgorithm.h>
#include <itkImageToImageFilter.h>
#include <itkIOCommon.h>
#include <itkVector.h>
#include <itkListSample.h>
#include <vnl/vnl_vector.h>
#include <itkFixedArray.h>
#include <itkMetaDataObject.h>
#include "itkFastMarchingCostFunction.h"

#include "itkDtiFastMarchingTrackingFilter.h"
#include <itkRegularStepGradientDescentOptimizer.h>
#include "GtractTypes.h"
#include <map>
#include <string>

#include <itkIndex.h>
#include <itkMath.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkPointSet.h>
#include <itkBlobSpatialObject.h>
#include <itkSceneSpatialObject.h>

#include <iostream>

#include <itkNumericTraits.h>
#include <algorithm>

namespace itk
{
/*
 *
 */
template <typename TTensorImageType, typename TAnisotropyImageType, typename TCostImageType, typename TMaskImageType>
DtiFastMarchingTrackingFilter<TTensorImageType, TAnisotropyImageType, TCostImageType, TMaskImageType>::
  DtiFastMarchingTrackingFilter()
{
  m_CostIP = CostIPType::New();
  m_CostFN = CostFunctionType::New();
  m_GradientOP = OptimizerType::New();
  m_StartPointThreshold = 0.3;
  m_AnisotropyThreshold = 0.2;
  m_MaxStepSize = 1.0;
  m_MinStepSize = 0.01;
  m_CosFunctionStepSize = 1.0;
  m_NumberOfIterations = 200;
  m_StartPoints.clear();
}

/*
 *
 */
template <typename TTensorImageType, typename TAnisotropyImageType, typename TCostImageType, typename TMaskImageType>
void
DtiFastMarchingTrackingFilter<TTensorImageType, TAnisotropyImageType, TCostImageType, TMaskImageType>::PrintSelf(
  std::ostream & os,
  Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Input Cost Image: " << m_CostImage.GetPointer() << std::endl;
  os << indent << "Input Anisotropy Image: " << m_AnisotropyImage.GetPointer() << std::endl;
}

template <typename TTensorImageType, typename TAnisotropyImageType, typename TCostImageType, typename TMaskImageType>
void
DtiFastMarchingTrackingFilter<TTensorImageType, TAnisotropyImageType, TCostImageType, TMaskImageType>::InitializeSeeds()
{
  // ////////////////////////////////////////////////////////////////////////
  // Initialize the seed points
  // ////////////////////////////////////////////////////////////////////////

  using ConstMaskIteratorType = itk::ImageRegionConstIterator<MaskImageType>;
  ConstMaskIteratorType maskIt(m_StartingRegion, m_StartingRegion->GetLargestPossibleRegion());
  int                   count = 0;
  for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt)
  {
    typename Self::ContinuousIndexType        seed;
    typename ConstMaskIteratorType::IndexType pos = maskIt.GetIndex();
    seed[0] = pos[0];
    seed[1] = pos[1];
    seed[2] = pos[2];
    const float ai = m_ScalarIP->EvaluateAtContinuousIndex(seed);
    if (maskIt.Get() && ai >= m_SeedThreshold)
    {
      m_StartPoints.push_back(startIndex);
      count++;
    }
  }
  std::cerr << "Number of Seeds: " << count << std::endl;
}

/*
 *
 */
template <typename TTensorImageType, typename TAnisotropyImageType, typename TCostImageType, typename TMaskImageType>
void
DtiFastMarchingTrackingFilter<TTensorImageType, TAnisotropyImageType, TCostImageType, TMaskImageType>::Update()
{
  this->m_Output = vtkPolyData::New();
  this->m_StartPoints.clear();
  this->m_ScalarIP->SetInputImage(this->m_AnisotropyImage);
  this->m_VectorIP->SetInputImage(this->m_TensorImage);
  this->m_EndIP->SetInputImage(this->m_EndingRegion);
  this->m_CostIP->SetInputImage(this->m_CostImage);

  this->InitializeSeeds();

  while (!m_StartPoints.empty())
  {
    Self::ContinuousIndexType inputIndex = m_StartPoints.front();
    m_StartPoints.pop_front();

    // check if index is within the input image region
    if (!costImage->GetBufferedRegion().IsInsideInWorldSpace(inputIndex))
    {
      continue;
    }

    // Get fiber through Gradient Descent
    this->GradientDescent(inputIndex);
  }
}

/*
 *
 */

template <typename TTensorImageType, typename TAnisotropyImageType, typename TCostImageType, typename TMaskImageType>
void
DtiFastMarchingTrackingFilter<TTensorImageType, TAnisotropyImageType, TCostImageType, TMaskImageType>::GradientDescent(
  CostIndexType & index)
{
  CostIndexType inputIndex = index;
  CostIndexType tmpIndex = index;
  FiberType     fiber; // list of fiberpoints

  fiber.clear();
  CostFiberPointType fp;
  float              anisotropy;
  // float ci;  //cost value
  bool   completeFiber = true;
  bool   pass = false;
  double gradientTol = 0.001;

  vtkPoints *     fiber = vtkPoints::New();
  vtkFloatArray * fiberTensors = vtkFloatArray::New();
  fiberTensors->SetName("Tensors");
  fiberTensors->SetNumberOfComponents(9);
  vtkFloatArray * fiberAnisotropy = vtkFloatArray::New();
  fiberAnisotropy->SetName("Anisotropy");
  vtkFloatArray * fiberAnisotropySum = vtkFloatArray::New();
  fiberAnisotropySum->SetName("Anisotropy-Sum");
  vtkFloatArray * fiberCost = vtkFloatArray::New();
  fiberCost->SetName("Cost");

  /*Set up the Cost function*/
  m_CostFN->SetCostImage(m_CostImage);

  /*Set up the gradient descent optimizer*/
  m_GradientOP->SetCostFunction(m_CostFN);

  unsigned int spaceDimension = m_CostFN->GetNumberOfParameters();

  ParametersType initialPosition(spaceDimension);
  ScalesType     parametersScale(spaceDimension);
  DerivativeType gradient(spaceDimension);
  for (unsigned int i = 0; i < spaceDimension; i++)
  {
    initialPosition[i] = inputIndex[i];
    parametersScale[i] = 1.0; // 1 by default
  }

  /* Set up rest of Optimizer parameters except intialPosition*/
  m_GradientOP->MaximizeOn();
  m_GradientOP->SetScales(parametersScale);
  m_GradientOP->SetMaximumStepLength(m_MaxStepSize);
  m_GradientOP->SetGradientMagnitudeTolerance(gradientTol);
  m_GradientOP->SetMinimumStepLength(m_MinStepSize);
  m_GradientOP->SetNumberOfIterations(1);
  m_GradientOP->SetRelaxationFactor(.8);

  /* Initial points are StartPoints and valid threshold has been done in InitializeStartPoints()
     Thus, we can add intitial points immediately to fiber*/
  anisotropy = this->m_ScalarIP->EvaluateAtContinuousIndex(inputIndex);
  anisotropySum = anisotropy;
  fiberAnisotropy->InsertNextValue(anisotropy);
  fiberAnisotropySum->InsertNextValue(anisotropySum);
  fiberCost->InsertNextValue(m_CostFN->GetValue(initialPosition));

  typename Self::PointType p;
  this->ContinuousIndexToMM(inputIndex, p);
  fiber->InsertNextPoint(p.GetDataPointer());

  typename Self::TensorImagePixelType tensorPixel = this->m_VectorIP->EvaluateAtContinuousIndex(index);

  TMatrix fullTensorPixel(3, 3);
  fullTensorPixel = Tensor2Matrix(tensorPixel);
  fiberTensors->InsertNextTupleValue(fullTensorPixel.data_block());

  m_GradientOP->SetInitialPosition(initialPosition);

  /* Start Optimization : Must Start Optimization before Resume Operation*/
  try
  {
    m_GradientOP->StartOptimization();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cout << "Exception thrown ! " << std::endl;
    std::cout << "An error ocurred during Optimization" << std::endl;
    std::cout << "Location    = " << e.GetLocation() << std::endl;
    std::cout << "Description = " << e.GetDescription() << std::endl;
  }

  ParametersType currentPosition(spaceDimension);
  unsigned int   count = 0;
  /*Resume Optimization */
  for (unsigned int j = 0; j < (m_NumberOfIterations - 1); j++)
  {
    pass = false;
    currentPosition = m_GradientOP->GetCurrentPosition();
    for (unsigned int k = 0; k < spaceDimension; k++) //
    {
      double diff = itk::Math::abs(currentPosition[k] - tmpIndex[k]);

      // if last point is repeated, stop Gradient Descent
      if (diff > gradientTol)
      {
        pass = true;
        tmpIndex[k] = currentPosition[k];
      }
    }

    anisotropy = m_AnisoIP->EvaluateAtContinuousIndex(tmpIndex);
    anisotropySum += anisotropy;

    if (pass)
    {
      if (anisotropy >= m_AnisotropyThreshold)
      {
        // Add current point (fiber point) to fiber
        this->ContinuousIndexToMM(tmpIndex, p);
        fiber->InsertNextPoint(p.GetDataPointer());
        fiberAnisotropy->InsertNextValue(anisotropy);
        fiberAnisotropySum->InsertNextValue(anisotropySum);
        fiberCost->InsertNextValue(m_CostFN->GetValue(currentPosition));

        // Add the Tensor to the Scalar Data
        tensorPixel = this->m_VectorIP->EvaluateAtContinuousIndex(index);
        fullTensorPixel = Tensor2Matrix(tensorPixel);
        fiberTensors->InsertNextTupleValue(fullTensorPixel.data_block());

        // Reset gradient optimizer with current point as starting point
        m_GradientOP->SetInitialPosition(currentPosition);
        m_GradientOP->SetNumberOfIterations(count + 2);

        // Do next iteration
        try
        {
          m_GradientOP->ResumeOptimization();
        }
        catch (itk::ExceptionObject & e)
        {
          std::cout << "Exception thrown ! " << std::endl;
          std::cout << "An error ocurred during Optimization" << std::endl;
          std::cout << "Location    = " << e.GetLocation() << std::endl;
          std::cout << "Description = " << e.GetDescription() << std::endl;
        }
      }
    }

    count += 1;
  } // end outer for loop

  gradient = m_GradientOP->GetGradient();
  for (unsigned int i = 0; i < spaceDimension; i++)
  {
    if ((gradient[i] > 0.0) || (gradient[i] < 0.0))
    {
      completeFiber = false;
    }
  }

  if (completeFiber)
  {
    this->AddFiberToOutput(fiber, fiberTensors);
  }
} // end Gradient Descent
} // namespace itk

#endif
