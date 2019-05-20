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
/*
 *  commandIterationupdate.h
 *  ThirionCLP
 *
 *  Created by Yong Qiang Zhao on 8/21/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include "itkCommand.h"
#include "itkMultiResolutionPDEDeformableRegistration.h"
#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include "itkDemonsRegistrationFilter.h"
#include "itkDisplacementFieldJacobianDeterminantFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkWarpHarmonicEnergyCalculator.h"
#include "itkVectorCentralDifferenceImageFunction.h"
#include "itkPDEDeformableRegistrationFilter.h"

#include "BRAINSCommonLib.h"
#include "DebugImageViewerClient.h"
#include "DebugImageWrite.h"
#include "GenericTransformImage.h"

#include "BRAINSDemonWarpTemplates.h"
#include "itkFastSymmetricForcesDemonsRegistrationFilter.h"
#include "itkDiffeomorphicDemonsRegistrationWithMaskFilter.h"

template < typename TPixel = float, unsigned int VImageDimension = 3 >
class CommandIterationUpdate : public itk::Command
{
public:
  using Self = CommandIterationUpdate;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer< Self >;

  using InternalImageType = itk::Image< TPixel, VImageDimension >;
  using VectorPixelType = itk::Vector< TPixel, VImageDimension >;
  using DisplacementFieldType = itk::Image< VectorPixelType, VImageDimension >;
  using VelocityFieldType = itk::Image< VectorPixelType, VImageDimension >;

  using DemonsRegistrationFilterType =
    itk::DemonsRegistrationFilter< InternalImageType, InternalImageType, DisplacementFieldType >;

  using DiffeomorphicDemonsRegistrationFilterType =
    itk::DiffeomorphicDemonsRegistrationFilter< InternalImageType, InternalImageType, DisplacementFieldType >;

  using DiffeomorphicDemonsRegistrationWithMaskFilterType =
    itk::DiffeomorphicDemonsRegistrationWithMaskFilter< InternalImageType, InternalImageType, DisplacementFieldType >;

  using FastSymmetricForcesDemonsRegistrationFilterType =
    itk::FastSymmetricForcesDemonsRegistrationFilter< InternalImageType, InternalImageType, DisplacementFieldType >;

  using MultiResRegistrationFilterType =
    itk::MultiResolutionPDEDeformableRegistration< InternalImageType, InternalImageType, DisplacementFieldType,
                                                   TPixel >;

  using JacobianFilterType =
    itk::DisplacementFieldJacobianDeterminantFilter< DisplacementFieldType, TPixel, InternalImageType >;

  using MinMaxFilterType = itk::MinimumMaximumImageCalculator< InternalImageType >;

  using HarmonicEnergyCalculatorType = itk::WarpHarmonicEnergyCalculator< DisplacementFieldType >;

  using WarpGradientCalculatorType = itk::VectorCentralDifferenceImageFunction< DisplacementFieldType >;

  using WarpGradientType = typename WarpGradientCalculatorType::OutputType;

  itkNewMacro( Self );

private:
  std::ofstream                                  m_Fid;
  bool                                           m_headerwritten;
  typename JacobianFilterType::Pointer           m_JacobianFilter;
  typename MinMaxFilterType::Pointer             m_Minmaxfilter;
  typename HarmonicEnergyCalculatorType::Pointer m_HarmonicEnergyCalculator;
  typename DisplacementFieldType::ConstPointer   m_TrueField;
  typename WarpGradientCalculatorType::Pointer   m_TrueWarpGradientCalculator;
  typename WarpGradientCalculatorType::Pointer   m_CompWarpGradientCalculator;
  typename InternalImageType::Pointer            m_MovingImage;
  typename InternalImageType::Pointer            m_FixedImage;

public:
  void
  SetMovingImage( typename InternalImageType::Pointer & img )
  {
    m_MovingImage = img;
  }

  void
  SetFixedImage( typename InternalImageType::Pointer & img )
  {
    m_FixedImage = img;
  }

  void
  SetTrueField( const DisplacementFieldType * truefield )
  {
    m_TrueField = truefield;

    m_TrueWarpGradientCalculator = WarpGradientCalculatorType::New();
    m_TrueWarpGradientCalculator->SetInputImage( m_TrueField );

    m_CompWarpGradientCalculator = WarpGradientCalculatorType::New();
  }

  void
  Execute( itk::Object * caller, const itk::EventObject & event ) override
  {
    Execute( (const itk::Object *)caller, event );
  }

  void
  Execute( const itk::Object * object, const itk::EventObject & event ) override
  {
    if ( !( itk::IterationEvent().CheckEvent( &event ) ) )
    {
      return;
    }

    typename DisplacementFieldType::Pointer deffield = nullptr;
    unsigned int                            iter = std::numeric_limits< unsigned int >::max();
    double                                  metricbefore = -1.0;

    if ( const DiffeomorphicDemonsRegistrationFilterType * DDfilter =
           dynamic_cast< const DiffeomorphicDemonsRegistrationFilterType * >( object ) )
    {
      iter = DDfilter->GetElapsedIterations() - 1;
      metricbefore = DDfilter->GetMetric();
      deffield = const_cast< DiffeomorphicDemonsRegistrationFilterType * >( DDfilter )->GetDisplacementField();
    }
    else if ( const DiffeomorphicDemonsRegistrationWithMaskFilterType * DDWMfilter =
                dynamic_cast< const DiffeomorphicDemonsRegistrationWithMaskFilterType * >( object ) )
    {
      iter = DDWMfilter->GetElapsedIterations() - 1;
      metricbefore = DDWMfilter->GetMetric();
      deffield =
        const_cast< DiffeomorphicDemonsRegistrationWithMaskFilterType * >( DDWMfilter )->GetDisplacementField();
    }
    else if ( const FastSymmetricForcesDemonsRegistrationFilterType * FSDfilter =
                dynamic_cast< const FastSymmetricForcesDemonsRegistrationFilterType * >( object ) )
    {
      iter = FSDfilter->GetElapsedIterations() - 1;
      metricbefore = FSDfilter->GetMetric();
      deffield = const_cast< FastSymmetricForcesDemonsRegistrationFilterType * >( FSDfilter )->GetDisplacementField();
    }
    else if ( const DemonsRegistrationFilterType * Dfilter =
                dynamic_cast< const DemonsRegistrationFilterType * >( object ) )
    {
      iter = Dfilter->GetElapsedIterations() - 1;
      metricbefore = Dfilter->GetMetric();
      deffield = const_cast< DemonsRegistrationFilterType * >( Dfilter )->GetDisplacementField();
    }
    else if ( const MultiResRegistrationFilterType * multiresfilter =
                dynamic_cast< const MultiResRegistrationFilterType * >( object ) )
    {
      std::cout << "Finished Multi-resolution iteration :" << multiresfilter->GetCurrentLevel() - 1 << std::endl;
      std::cout << "==============================" << std::endl << std::endl;
    }
    else
    {
      return;
    }

    if ( deffield.IsNotNull() )
    {
      std::cout << iter << ": MSE " << metricbefore << " - ";

      double fieldDist = -1.0;
      double fieldGradDist = -1.0;
      double tmp;
      if ( m_TrueField )
      {
        using FieldIteratorType = itk::ImageRegionConstIteratorWithIndex< DisplacementFieldType >;
        FieldIteratorType currIter( deffield, deffield->GetLargestPossibleRegion() );
        FieldIteratorType trueIter( m_TrueField, deffield->GetLargestPossibleRegion() );

        m_CompWarpGradientCalculator->SetInputImage( deffield );

        fieldDist = 0.0;
        fieldGradDist = 0.0;
        for ( currIter.GoToBegin(), trueIter.GoToBegin(); !currIter.IsAtEnd(); ++currIter, ++trueIter )
        {
          fieldDist += ( currIter.Value() - trueIter.Value() ).GetSquaredNorm();

          // No need to add Id matrix here as we do a substraction
          tmp = ( ( m_CompWarpGradientCalculator->EvaluateAtIndex( currIter.GetIndex() ) -
                    m_TrueWarpGradientCalculator->EvaluateAtIndex( trueIter.GetIndex() ) )
                    .GetVnlMatrix() )
                  .frobenius_norm();
          fieldGradDist += tmp * tmp;
        }
        fieldDist = sqrt( fieldDist / (double)( deffield->GetLargestPossibleRegion().GetNumberOfPixels() ) );
        fieldGradDist = sqrt( fieldGradDist / (double)( deffield->GetLargestPossibleRegion().GetNumberOfPixels() ) );

        std::cout << "d(.,true) " << fieldDist << " - ";
        std::cout << "d(.,Jac(true)) " << fieldGradDist << " - ";
      }
#if defined( USE_DebugImageViewer )
      if ( DebugImageDisplaySender.Enabled() )
      {
        DebugImageDisplaySender.SendImage< DisplacementFieldType >( deffield, 0, 0 );
        DebugImageDisplaySender.SendImage< DisplacementFieldType >( deffield, 1, 1 );
        DebugImageDisplaySender.SendImage< DisplacementFieldType >( deffield, 2, 2 );
        typename InternalImageType::Pointer DeformedMovingImagePtr =
          TransformWarp< InternalImageType, InternalImageType, DisplacementFieldType >(
            m_MovingImage, deffield, 0, GetInterpolatorFromString< InternalImageType >( "Linear" ), deffield );
        DebugImageDisplaySender.SendImage< InternalImageType >( DeformedMovingImagePtr, 3 );
        //        std::cerr << std::endl << "************IMAGES
        // SENT*************" << std::endl;
      }
#endif // defined(USE_DebugImageViewer)

      m_HarmonicEnergyCalculator->SetImage( deffield );
      m_HarmonicEnergyCalculator->Compute();
      const double harmonicEnergy = m_HarmonicEnergyCalculator->GetHarmonicEnergy();
      std::cout << "harmo. " << harmonicEnergy << " - ";

      m_JacobianFilter->SetInput( deffield );
      m_JacobianFilter->UpdateLargestPossibleRegion();

      const unsigned int numPix = m_JacobianFilter->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();

      TPixel * pix_start = m_JacobianFilter->GetOutput()->GetBufferPointer();
      TPixel * pix_end = pix_start + numPix;

      TPixel * jac_ptr;

      // Get percentage of det(Jac) below 0
      unsigned int jacBelowZero( 0u );
      for ( jac_ptr = pix_start; jac_ptr != pix_end; ++jac_ptr )
      {
        if ( *jac_ptr <= 0.0 )
        {
          ++jacBelowZero;
        }
      }
      const double jacBelowZeroPrc = static_cast< double >( jacBelowZero ) / static_cast< double >( numPix );

      // Get min an max jac
      const double minJac = *( std::min_element( pix_start, pix_end ) );
      const double maxJac = *( std::max_element( pix_start, pix_end ) );

      // Get some quantiles
      // We don't need the jacobian image
      // we can modify/sort it in place
      jac_ptr = pix_start + static_cast< unsigned int >( 0.002 * numPix );
      std::nth_element( pix_start, jac_ptr, pix_end );
      const double Q002 = *jac_ptr;

      jac_ptr = pix_start + static_cast< unsigned int >( 0.01 * numPix );
      std::nth_element( pix_start, jac_ptr, pix_end );
      const double Q01 = *jac_ptr;

      jac_ptr = pix_start + static_cast< unsigned int >( 0.99 * numPix );
      std::nth_element( pix_start, jac_ptr, pix_end );
      const double Q99 = *jac_ptr;

      jac_ptr = pix_start + static_cast< unsigned int >( 0.998 * numPix );
      std::nth_element( pix_start, jac_ptr, pix_end );
      const double Q998 = *jac_ptr;

      std::cout << "max|Jac| " << maxJac << " - "
                << "min|Jac| " << minJac << " - "
                << "ratio(|Jac|<=0) " << jacBelowZeroPrc << std::endl;

      if ( this->m_Fid.is_open() )
      {
        if ( !m_headerwritten )
        {
          this->m_Fid << "Iteration"
                      << ", MSE before"
                      << ", Harmonic energy"
                      << ", min|Jac|"
                      << ", 0.2% |Jac|"
                      << ", 01% |Jac|"
                      << ", 99% |Jac|"
                      << ", 99.8% |Jac|"
                      << ", max|Jac|"
                      << ", ratio(|Jac|<=0)";

          if ( m_TrueField )
          {
            this->m_Fid << ", dist(warp,true warp)"
                        << ", dist(Jac,true Jac)";
          }

          this->m_Fid << std::endl;

          m_headerwritten = true;
        }

        this->m_Fid << iter << ", " << metricbefore << ", " << harmonicEnergy << ", " << minJac << ", " << Q002 << ", "
                    << Q01 << ", " << Q99 << ", " << Q998 << ", " << maxJac << ", " << jacBelowZeroPrc;

        if ( m_TrueField )
        {
          this->m_Fid << ", " << fieldDist << ", " << fieldGradDist;
        }

        this->m_Fid << std::endl;
      }
    }
  }

protected:
  CommandIterationUpdate()
    : m_Fid( "metricvalues.csv" )
    , m_headerwritten( false )
  {
    m_JacobianFilter = JacobianFilterType::New();
    m_JacobianFilter->SetUseImageSpacing( true );
    m_JacobianFilter->ReleaseDataFlagOn();

    m_Minmaxfilter = MinMaxFilterType::New();

    m_HarmonicEnergyCalculator = HarmonicEnergyCalculatorType::New();

    m_TrueField = nullptr;
    m_TrueWarpGradientCalculator = nullptr;
    m_CompWarpGradientCalculator = nullptr;
  }

  ~CommandIterationUpdate() override { this->m_Fid.close(); }
};
