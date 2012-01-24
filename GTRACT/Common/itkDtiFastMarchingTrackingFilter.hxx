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

#ifndef __itkDtiFastMarchingTrackingFilter_hxx
#define __itkDtiFastMarchingTrackingFilter_hxx

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
#include <vnl/vnl_math.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkPointSet.h>

#include <iostream>

#include <itkNumericTraits.h>
#include <algorithm>

namespace itk
{
/*
 *
 */
template <class TTensorImageType, class TAnisotropyImageType, class TCostImageType, class TMaskImageType>
DtiFastMarchingTrackingFilter<TTensorImageType, TAnisotropyImageType, TCostImageType, TMaskImageType>
::DtiFastMarchingTrackingFilter()
{
  m_CostIP     = CostIPType::New();
  m_CostFN     = CostFunctionType::New();
  m_GradientOP = OptimizerType::New();
  m_MaxStepSize          = 1.0;
  m_MinStepSize          = 0.01;
  m_CostFunctionStepSize = 1.0;
  m_NumberOfIterations   = 200;
  m_StartPoints.clear();
}

/*
 *
 */
template <class TTensorImageType, class TAnisotropyImageType, class TCostImageType, class TMaskImageType>
void
DtiFastMarchingTrackingFilter<TTensorImageType, TAnisotropyImageType, TCostImageType, TMaskImageType>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Max Step Size: " << m_MaxStepSize  << std::endl;
  os << indent << "Min Step Size: " << m_MinStepSize  << std::endl;
  os << indent << "Cost Function Step Size: " << m_CostFunctionStepSize  << std::endl;
  os << indent << "Number of Iterations: " << m_NumberOfIterations  << std::endl;
}

template <class TTensorImageType, class TAnisotropyImageType, class TCostImageType, class TMaskImageType>
void
DtiFastMarchingTrackingFilter<TTensorImageType, TAnisotropyImageType, TCostImageType, TMaskImageType>
::InitializeSeeds()
{
  // ////////////////////////////////////////////////////////////////////////
  // Initialize the seed points
  // ////////////////////////////////////////////////////////////////////////
  typedef itk::ImageRegionConstIterator<TMaskImageType> ConstMaskIteratorType;
  ConstMaskIteratorType maskIt( this->m_StartingRegion, this->m_StartingRegion->GetLargestPossibleRegion() );

  int count = 0;
  for( maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt )
    {
    typename Self::ContinuousIndexType        seed;
    typename ConstMaskIteratorType::IndexType pos = maskIt.GetIndex();
    seed[0] = pos[0];  seed[1] = pos[1];  seed[2] = pos[2];
    const float ai = this->m_ScalarIP->EvaluateAtContinuousIndex(seed);
    if( maskIt.Get()  &&  ai >= this->m_SeedThreshold )
      {
      m_StartPoints.push_back( seed );
      count++;
      }
    }
  std::cerr << "Number of Seeds: " << count << std::endl;
}

/*
 *
 */
template <class TTensorImageType, class TAnisotropyImageType, class TCostImageType, class TMaskImageType>
void
DtiFastMarchingTrackingFilter<TTensorImageType, TAnisotropyImageType, TCostImageType, TMaskImageType>
::Update()
{
  this->m_Output = vtkPolyData::New();
  this->m_StartPoints.clear();
  this->m_ScalarIP->SetInputImage(this->m_AnisotropyImage);
  this->m_VectorIP->SetInputImage(this->m_TensorImage);
  this->m_EndIP->SetInputImage(this->m_EndingRegion);
  this->m_CostIP->SetInputImage(this->m_CostImage);

  this->InitializeSeeds();

  while( !m_StartPoints.empty() )
    {
    typename Self::ContinuousIndexType inputIndex = m_StartPoints.front();
    m_StartPoints.pop_front();
    // std::cout << "Tracking point: " << inputIndex << std::endl;
    // check if index is within the input image region
    if( !m_CostImage->GetBufferedRegion().IsInside( inputIndex ) )
      {
      continue;
      }

    // Get fiber through Gradient Descent
    this->GradientDescent( inputIndex );
    }
}

/*
 *
 */

template <class TTensorImageType, class TAnisotropyImageType, class TCostImageType, class TMaskImageType>
void
DtiFastMarchingTrackingFilter<TTensorImageType, TAnisotropyImageType, TCostImageType, TMaskImageType>
::GradientDescent( ContinuousIndexType & index)
{
  typename Self::ContinuousIndexType inputIndex = index;
  typename Self::ContinuousIndexType tmpIndex = index;
  float  anisotropy;
  float  anisotropySum = 0;
  bool   completeFiber = true;
  bool   pass = false;
  double gradientTol = 0.001;

  vtkPoints *    fiber = vtkPoints::New();
  vtkFloatArray *fiberTensors = vtkFloatArray::New();

  fiberTensors->SetName("Tensors");
  fiberTensors->SetNumberOfComponents(9);
  vtkFloatArray *fiberAnisotropy = vtkFloatArray::New();
  fiberAnisotropy->SetName("Anisotropy");
  vtkFloatArray *fiberAnisotropySum = vtkFloatArray::New();
  fiberAnisotropySum->SetName("Anisotropy-Sum");
  vtkFloatArray *fiberCost = vtkFloatArray::New();
  fiberCost->SetName("Cost");

  /*Set up the Cost function*/
  m_CostFN->SetCostImage( m_CostImage );

  /*Set up the gradient descent optimizer*/
  m_GradientOP->SetCostFunction( m_CostFN );

  unsigned int spaceDimension = m_CostFN->GetNumberOfParameters();

  ParametersType initialPosition( spaceDimension );

  ScalesType parametersScale( spaceDimension );

  DerivativeType gradient( spaceDimension );
  for( unsigned int i = 0; i < spaceDimension; i++ )
    {
    initialPosition[i] = inputIndex[i];
    parametersScale[i] = 1.0; // 1 by default
    }

  /* Set up rest of Optimizer parameters except intialPosition*/
  m_GradientOP->MaximizeOn();
  m_GradientOP->SetScales( parametersScale );
  m_GradientOP->SetMaximumStepLength(m_MaxStepSize);
  m_GradientOP->SetGradientMagnitudeTolerance( gradientTol );
  m_GradientOP->SetMinimumStepLength( m_MinStepSize );
  m_GradientOP->SetNumberOfIterations( 1 );
  m_GradientOP->SetRelaxationFactor( .8 );

  /* Initial points are StartPoints and valid threshold has been done in InitializeStartPoints()
     Thus, we can add intitial points immediately to fiber*/
  anisotropy = this->m_ScalarIP->EvaluateAtContinuousIndex( inputIndex );
  anisotropySum = anisotropy;
  fiberAnisotropy->InsertNextValue( anisotropy );
  fiberAnisotropySum->InsertNextValue( anisotropySum );
  fiberCost->InsertNextValue( m_CostFN->GetValue(initialPosition) );

  typename Self::PointType p;
  this->ContinuousIndexToMM(inputIndex, p);
  fiber->InsertNextPoint( p.GetDataPointer() );

  typename Self::TensorImagePixelType tensorPixel = this->m_VectorIP->EvaluateAtContinuousIndex(index);

  TMatrix fullTensorPixel(3, 3); fullTensorPixel = Tensor2Matrix( tensorPixel );
  fiberTensors->InsertNextTupleValue( fullTensorPixel.data_block() );

  m_GradientOP->SetInitialPosition( initialPosition );

  /* Start Optimization : Must Start Optimization before Resume Operation*/
  try
    {
    m_GradientOP->StartOptimization();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception thrown ! " << std::endl;
    std::cout << "An error ocurred during Optimization" << std::endl;
    std::cout << "Location    = " << e.GetLocation()    << std::endl;
    std::cout << "Description = " << e.GetDescription() << std::endl;
    }

  ParametersType currentPosition( spaceDimension );

  unsigned int count = 0;
  /*Resume Optimization */
  for( unsigned int j = 0; j < ( m_NumberOfIterations - 1 ); j++ )
    {
    pass = false;
    currentPosition = m_GradientOP->GetCurrentPosition();
    for( unsigned int k = 0; k < spaceDimension; k++ )  //
      {
      double diff = vnl_math_abs(currentPosition[k] - tmpIndex[k]);

      // if last point is repeated, stop Gradient Descent
      if( diff > gradientTol )
        {
        pass = true;
        tmpIndex[k] = currentPosition[k];
        }
      }

    anisotropy = this->m_ScalarIP->EvaluateAtContinuousIndex(tmpIndex);
    anisotropySum += anisotropy;

    if( pass )
      {
      if( anisotropy >= this->m_AnisotropyThreshold )
        {
        // Add current point (fiber point) to fiber
        this->ContinuousIndexToMM(tmpIndex, p);
        fiber->InsertNextPoint( p.GetDataPointer() );
        fiberAnisotropy->InsertNextValue( anisotropy );
        fiberAnisotropySum->InsertNextValue( anisotropySum );
        fiberCost->InsertNextValue( m_CostFN->GetValue(currentPosition) );

        // Add the Tensor to the Scalar Data
        tensorPixel = this->m_VectorIP->EvaluateAtContinuousIndex(index);
        fullTensorPixel = Tensor2Matrix( tensorPixel );
        fiberTensors->InsertNextTupleValue( fullTensorPixel.data_block() );

        // Reset gradient optimizer with current point as starting point
        m_GradientOP->SetInitialPosition( currentPosition );
        m_GradientOP->SetNumberOfIterations(count + 2);

        // Do next iteration
        try
          {
          m_GradientOP->ResumeOptimization();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception thrown ! " << std::endl;
          std::cout << "An error ocurred during Optimization" << std::endl;
          std::cout << "Location    = " << e.GetLocation()    << std::endl;
          std::cout << "Description = " << e.GetDescription() << std::endl;
          }
        }
      else
        {
        // std::cout << "Terminate Low FA" << std::endl;
        }
      }

    count += 1;
    } // end outer for loop

  gradient = m_GradientOP->GetGradient();
  for( unsigned int i = 0; i < spaceDimension; i++ )
    {
    if( ( gradient[i] > 0.0 ) || ( gradient[i] < 0.0 ) )
      {
      completeFiber = false;
      }
    }
  // std::cout << "Gradient: " << gradient[0] << " " << gradient[1] << " " <<
  // gradient[2] << std::endl;
  // std::cout << "Fiber COmplete: " << completeFiber << std::endl;
  if( completeFiber )
    {
    this->AddFiberToOutput( fiber, fiberTensors );
    }
} // end Gradient Descent
}   // namespace itk

#endif
