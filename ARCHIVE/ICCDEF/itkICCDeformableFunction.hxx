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
 *  itkICCDeformableFunction.txx
 *  ICCDeformationTools
 *
 *  Created by Yongqiang Zhao on 4/21/09.
 *  Copyright 2009 The University of Iowa. All rights reserved.
 *
 */

#include "itkICCDeformableFunction.h"

#ifndef __itkICCDeformableFunction_txx
#  define __itkICCDeformableFunction_txx

#  include <sstream>

#  include "itkICCDeformableFunction.h"
#  include "itkExceptionObject.h"
#  include "itkMath.h"
// #include "itkIterativeInverseDisplacementFieldImageFilter.h"
// #include "itkIterativeInverseDisplacementFieldImageFilter1.h"
#  include "itkICCIterativeInverseDisplacementFieldImageFilter.h"
#  include "itkSubtractImageFilter.h"
#  include "itkAddImageFilter.h"
#  include "itkImageFileWriter.h"
#  include "itkImageMaskSpatialObject.h"
#  include "itkVectorIndexSelectionCastImageFilter.h"
#  include "itkGaussianOperator.h"
#  include "itkVectorNeighborhoodOperatorImageFilter.h"
#  include "itkConstNeighborhoodIterator.h"
#  include "itkMaskImageFilter.h"
#  include "itkLabelStatisticsImageFilter.h"
#  include "itkThinPlateSplineKernelTransform.h"

namespace itk
{
/**
 * Default constructor
 */
template < typename TFixedImage, typename TMovingImage, typename TDisplacementField >
ICCDeformableFunction< TFixedImage, TMovingImage, TDisplacementField >::ICCDeformableFunction()
{
  RadiusType   r;
  unsigned int j;

  for ( j = 0; j < ImageDimension; j++ )
  {
    r[j] = 0;
  }
  this->SetRadius( r );

  m_TimeStep = 1.0;
  m_DenominatorThreshold = 1e-9;
  m_IntensityDifferenceThreshold = 0.001;
  m_MaximumUpdateStepLength = 0.5;

  this->SetMovingImage( nullptr );
  this->SetFixedImage( nullptr );
  m_FixedImageSpacing.Fill( 1.0 );
  m_FixedImageOrigin.Fill( 0.0 );
  m_FixedImageDirection.SetIdentity();
  m_Normalizer = 0.0;

  this->m_UseGradientType = Symmetric;

  typename DefaultMovingInterpolatorType::Pointer interp = DefaultMovingInterpolatorType::New();

  m_MovingImageInterpolator = static_cast< InterpolatorMovingType * >( interp.GetPointer() );

  m_MovingImageWarper = WarperMovingType::New();
  m_MovingImageWarper->SetInterpolator( m_MovingImageInterpolator );
  m_MovingMaskImageWarper = MaskWarperType::New();

  //  m_MovingImageWarper->SetEdgePaddingValue( NumericTraits<MovingPixelType>::max() );

  m_Metric = NumericTraits< double >::max();
  m_SumOfSquaredDifference = 0.0;
  m_NumberOfPixelsProcessed = 0L;
  m_RMSChange = NumericTraits< double >::max();
  m_SumOfSquaredChange = 0.0;

  m_UpdateBuffer = DisplacementFieldType::New();
  m_InverseUpdateBuffer = DisplacementFieldType::New();
  m_Coefficient = DisplacementFieldFFTType::New();
  m_WarpedMaskImage = MaskImageType::New();

  m_SimilarityWeight = 1.0;
  m_LandmarkWeight = 0.0;
  m_FixedImageBackground = 70.0;
  m_MovingImageBackground = 70.0;
  m_BackgroundFilledValue = 35.0;
  m_FixedLandmark = PointSetType::New();
  m_MovingLandmark = PointSetType::New();
  m_UseConsistentIntensity = false;
  m_UseConsistentLandmark = false;
}

/*
 * Standard "PrintSelf" method.
 */
template < typename TFixedImage, typename TMovingImage, typename TDisplacementField >
void
ICCDeformableFunction< TFixedImage, TMovingImage, TDisplacementField >::PrintSelf( std::ostream & os,
                                                                                   Indent         indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "UseGradientType: ";
  os << m_UseGradientType << std::endl;
  os << indent << "MaximumUpdateStepLength: ";
  os << m_MaximumUpdateStepLength << std::endl;

  os << indent << "MovingImageIterpolator: ";
  os << m_MovingImageInterpolator.GetPointer() << std::endl;
  os << indent << "DenominatorThreshold: ";
  os << m_DenominatorThreshold << std::endl;
  os << indent << "IntensityDifferenceThreshold: ";
  os << m_IntensityDifferenceThreshold << std::endl;

  os << indent << "Metric: ";
  os << m_Metric << std::endl;
  os << indent << "SumOfSquaredDifference: ";
  os << m_SumOfSquaredDifference << std::endl;
  os << indent << "NumberOfPixelsProcessed: ";
  os << m_NumberOfPixelsProcessed << std::endl;
  os << indent << "RMSChange: ";
  os << m_RMSChange << std::endl;
  os << indent << "SumOfSquaredChange: ";
  os << m_SumOfSquaredChange << std::endl;
}

/**
 *
 */
template < typename TFixedImage, typename TMovingImage, typename TDisplacementField >
void
ICCDeformableFunction< TFixedImage, TMovingImage, TDisplacementField >::SetIntensityDifferenceThreshold(
  double threshold )
{
  m_IntensityDifferenceThreshold = threshold;
}

/**
 *
 */
template < typename TFixedImage, typename TMovingImage, typename TDisplacementField >
double
ICCDeformableFunction< TFixedImage, TMovingImage, TDisplacementField >::GetIntensityDifferenceThreshold() const
{
  return m_IntensityDifferenceThreshold;
}

/**
 * Set the function state values before each iteration
 */
template < typename TFixedImage, typename TMovingImage, typename TDisplacementField >
void
ICCDeformableFunction< TFixedImage, TMovingImage, TDisplacementField >::InitializeIteration()
{
  if ( !this->GetMovingImage() || !this->GetFixedImage() || !m_MovingImageInterpolator )
  {
    itkExceptionMacro( << "MovingImage, FixedImage and/or Interpolator not set" );
  }

  m_SumOfSquaredDifference = 0.0;
  m_NumberOfPixelsProcessed = 0L;
  m_SumOfSquaredChange = 0.0;

  // cache fixed image information
  m_FixedImageOrigin = this->GetFixedImage()->GetOrigin();
  m_FixedImageSpacing = this->GetFixedImage()->GetSpacing();
  m_FixedImageDirection = this->GetFixedImage()->GetDirection();
  // std::cout<<"Function!"<<std::endl;
  using ImageMaskSpatialObjectType = ImageMaskSpatialObject< Self::ImageDimension >;
  if ( this->GetMovingImageMask() && this->GetFixedImageMask() )
  {
    m_MovingMaskImageWarper->SetOutputOrigin( this->m_FixedImageOrigin );
    m_MovingMaskImageWarper->SetOutputSpacing( this->m_FixedImageSpacing );
    m_MovingMaskImageWarper->SetOutputDirection( this->m_FixedImageDirection );
    m_MovingMaskImageWarper->SetInput(
      dynamic_cast< ImageMaskSpatialObjectType * >( m_MovingMask.GetPointer() )->GetImage() );
    m_MovingMaskImageWarper->SetDisplacementField( this->GetDisplacementField() );
    m_MovingMaskImageWarper->GetOutput()->SetRequestedRegion( this->GetDisplacementField()->GetRequestedRegion() );
    m_MovingMaskImageWarper->Update();
    m_WarpedMaskImage = m_MovingMaskImageWarper->GetOutput();
    // itkUtil::WriteImage<MaskImageType>(m_WarpedMaskImage, "deformedMask.nii.gz");
  }

  using MaskFilterType = MaskImageFilter< MovingImageType, MaskImageType, MovingImageType >;
  typename MaskFilterType::Pointer mask = MaskFilterType::New();

  // Warped the moving and fixed image with two deformation fields
  m_MovingImageWarper->SetOutputOrigin( this->m_FixedImageOrigin );
  m_MovingImageWarper->SetOutputSpacing( this->m_FixedImageSpacing );
  m_MovingImageWarper->SetOutputDirection( this->m_FixedImageDirection );
  m_MovingImageWarper->SetInput( this->GetMovingImage() );
  m_MovingImageWarper->SetDisplacementField( this->GetDisplacementField() );
  m_MovingImageWarper->GetOutput()->SetRequestedRegion( this->GetDisplacementField()->GetRequestedRegion() );
  try
  {
    m_MovingImageWarper->Update();
  }
  catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
  }

  // Subtract fixed image from the Warped moving image
  typename SubtractImageType::Pointer subtract = SubtractImageType::New();

  if ( this->GetUseConsistentIntensity() )
  {
    typename DerivativeType::Pointer derivative = DerivativeType::New();

    derivative->SetInput( m_MovingImageWarper->GetOutput() );
    // compute the derivative of Warped image
    derivative->SetDirection( 0 );
    derivative->Update();

    typename FloatImageType::Pointer tempX = derivative->GetOutput();
    tempX->DisconnectPipeline();

    derivative->SetDirection( 1 );
    derivative->Update();

    typename FloatImageType::Pointer tempY = derivative->GetOutput();
    tempY->DisconnectPipeline();

    derivative->SetDirection( 2 );
    derivative->Update();

    typename FloatImageType::Pointer tempZ = derivative->GetOutput();
    tempZ->DisconnectPipeline();

    //  m_SimilarityWeight *= 4.0 * m_MaximumUpdateStepLength;
    //  typename TFixedImage::SizeType size = this->GetFixedImage()->GetLargestPossibleRegion().GetSize();

    //  std::cout<<"m_SimilarityWeight:"<<m_SimilarityWeight<<std::endl;
    //  m_SimilarityWeight = 0.1; //40.0/256.0;

    // Compute the similarity
    DisplacementFieldTypePointer similarity = DisplacementFieldType::New();
    similarity->SetRegions( this->GetDisplacementField()->GetLargestPossibleRegion() );
    similarity->Allocate();
    similarity->SetSpacing( this->GetDisplacementField()->GetSpacing() );
    similarity->SetOrigin( this->GetDisplacementField()->GetOrigin() );
    similarity->SetDirection( this->GetDisplacementField()->GetDirection() );

    IterationDisplacementFieldType def12Iter( similarity, similarity->GetRequestedRegion() );
    IterationImageType             sub12Iter( subtract->GetOutput(), subtract->GetOutput()->GetRequestedRegion() );
    ConstIterationImageType        fit( this->GetFixedImage(), this->GetFixedImage()->GetRequestedRegion() );
    IterationImageType mit( m_MovingImageWarper->GetOutput(), m_MovingImageWarper->GetOutput()->GetRequestedRegion() );
    IterationImageType xIter( tempX, tempX->GetRequestedRegion() );
    IterationImageType yIter( tempY, tempY->GetRequestedRegion() );
    IterationImageType zIter( tempZ, tempZ->GetRequestedRegion() );

    typename TDisplacementField::PixelType pixel;
    for ( def12Iter.GoToBegin(),
          fit.GoToBegin(),
          mit.GoToBegin(),
          xIter.GoToBegin(),
          yIter.GoToBegin(),
          zIter.GoToBegin();
          !def12Iter.IsAtEnd();
          ++def12Iter, ++fit, ++mit, ++xIter, ++yIter, ++zIter )
    {
      pixel[0] = 4.0 * m_MaximumUpdateStepLength * m_SimilarityWeight * ( mit.Get() - fit.Get() ) * xIter.Get();
      pixel[1] = 4.0 * m_MaximumUpdateStepLength * m_SimilarityWeight * ( mit.Get() - fit.Get() ) * yIter.Get();
      pixel[2] = 4.0 * m_MaximumUpdateStepLength * m_SimilarityWeight * ( mit.Get() - fit.Get() ) * zIter.Get();
      def12Iter.Set( pixel );
    }
    similarity->Modified();
    // Write the similarity measure image
    /*
      typename  DisplacementFieldFFTType::Pointer coeff =  DisplacementFieldFFTType::New();
       typename FFTWRealToComplexType::Pointer fft12 = FFTWRealToComplexType::New();
         fft12->SetInput(this->GetDisplacementField());
       fft12->Update();
       coeff=fft12->GetOutput();
       coeff->DisconnectPipeline();
    */

    typename FFTWRealToComplexType::Pointer fft_filter = FFTWRealToComplexType::New();
    fft_filter->SetInput( similarity );
    fft_filter->Update();

    ImageRegionIterator< DisplacementFieldFFTType > coeffsIter0( m_Coefficient, m_Coefficient->GetRequestedRegion() );
    ImageRegionConstIterator< ComplexImageType >    smoothIter( m_SmoothFilter, m_SmoothFilter->GetRequestedRegion() );
    ImageRegionConstIterator< DisplacementFieldFFTType > iter0( fft_filter->GetOutput(),
                                                                fft_filter->GetOutput()->GetRequestedRegion() );
    for ( iter0.GoToBegin(), smoothIter.GoToBegin(), coeffsIter0.GoToBegin(); !iter0.IsAtEnd();
          ++iter0, ++smoothIter, ++coeffsIter0 )
    {
      typename DisplacementFieldFFTType::PixelType pixel1 = coeffsIter0.Get();
      for ( unsigned i = 0; i < 3; i++ )
      {
        pixel1[i] = pixel1[i] - ( smoothIter.Get() * iter0.Get()[i] );
      }
      coeffsIter0.Set( pixel1 );
    }

    m_Coefficient->Modified();
  }

  if ( this->GetUseConsistentLandmark() )
  {
#  if 1
    // std::cout<<"landmarkWeight= "<<m_LandmarkWeight<<std::endl;
    if ( this->m_LandmarkWeight != 0 )
    {
      // Warp moving landmark: error
      // Kerneltransform needs pointset, thus, get points from landmark, then warp points, save point to pointset
      // first read points from Landmark
      // typename LandmarkType::PointListType list;
      using KernelTransformType = ThinPlateSplineKernelTransform< float, ImageDimension >;
      using KernelPointSetType = typename KernelTransformType::PointSetType;
      typename DisplacementFieldType::IndexType                   pixelIndex;
      typename DisplacementFieldType::IndexType                   pixelIndex1;
      typename MovingImageType::PointType                         movingPoint;
      typename MovingImageType::PointType                         warpedMovingPoint;
      typename MovingImageType::PointType                         fixedPoint;
      typename MovingImageType::PointType                         difference;
      typename KernelPointSetType::Pointer                        warpedMovingLandmark = KernelPointSetType::New();
      typename KernelPointSetType::Pointer                        fixedLandmark = KernelPointSetType::New();
      typename KernelPointSetType::PointsContainer::ConstIterator mlit;
      typename KernelPointSetType::PointsContainer::ConstIterator flit;

      // this->GetMovingImage()->Print(std::cout, 6);
      // this->GetFixedImage()->Print(std::cout,6);
      // std::cout<<"Number of Points:"<<m_MovingLandmark->GetNumberOfPoints()<<std::endl;
      float lmk_error_allowed = 0.0;
      mlit = m_MovingLandmark->GetPoints()->Begin();
      flit = m_FixedLandmark->GetPoints()->Begin();
      unsigned int pointID = 0;
      while ( mlit != m_MovingLandmark->GetPoints()->End() )
      {
        //    std::cout<<"Moving Landmark Position: "<<mlit.Value()<<std::endl;
        this->GetMovingImage()->TransformPhysicalPointToIndex( mlit.Value(), pixelIndex );
        this->GetFixedImage()->TransformIndexToPhysicalPoint( pixelIndex, movingPoint );
        //    std::cout<<"moving point: "<<movingPoint<<std::endl;
        this->GetFixedImage()->TransformPhysicalPointToIndex( flit.Value(), pixelIndex1 );
        this->GetFixedImage()->TransformIndexToPhysicalPoint( pixelIndex1, fixedPoint );
        //    std::cout<<"fixed point: "<<fixedPoint<<std::endl;
        for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
          warpedMovingPoint[j] = movingPoint[j] + this->GetDisplacementField()->GetPixel( pixelIndex )[j];
          difference[j] = fixedPoint[j] - warpedMovingPoint[j];
        }
#    if 1
        if ( std::abs( difference[0] ) > lmk_error_allowed || std::abs( difference[1] ) > lmk_error_allowed ||
             std::abs( difference[2] ) > lmk_error_allowed )
        {
          for ( unsigned int j = 0; j < ImageDimension; j++ )
          {
            warpedMovingPoint[j] = movingPoint[j] + difference[j];
          }
        }
#    endif
        this->GetFixedImage()->TransformPhysicalPointToIndex( warpedMovingPoint, pixelIndex1 );
        this->GetFixedImage()->TransformIndexToPhysicalPoint( pixelIndex1, warpedMovingPoint );
        warpedMovingLandmark->SetPoint( pointID, warpedMovingPoint );
        // std::cout<<"difference: "<<difference<<std::endl;
        fixedLandmark->SetPoint( pointID, movingPoint );
        ++mlit;
        ++flit;
        pointID++;
      }

      // warpedMovingLandmark->SetPoints(list);
      // Wondering to use physical space point or index point??????
      // Using kernel transform
      typename KernelTransformType::Pointer transform = KernelTransformType::New();
      transform->SetTargetLandmarks( fixedLandmark );
      transform->SetSourceLandmarks( warpedMovingLandmark );
      // transform->SetSourceLandmarks(fixedLandmark);
      // transform->SetTargetLandmarks(warpedMovingLandmark);
      transform->ComputeWMatrix();
      // std::cout<<"Displacement: "<<transform->GetDisplacements()<<std::endl;
      // ** Compute the displacement **//
      typename MovingImageType::PointType warpedPoint;
      typename MovingImageType::PointType point;
      IterationDisplacementFieldType      defIter( this->GetDisplacementField(),
                                              this->GetDisplacementField()->GetRequestedRegion() );
      ConstIterationImageType             it_m( this->GetMovingImage(), this->GetMovingImage()->GetRequestedRegion() );
      for ( it_m.GoToBegin(), defIter.GoToBegin(); !it_m.IsAtEnd(); ++it_m, ++defIter )
      {
        pixelIndex = it_m.GetIndex();
        this->GetFixedImage()->TransformIndexToPhysicalPoint( pixelIndex, point );
        warpedPoint = transform->TransformPoint( point );
        typename DisplacementFieldType::PixelType def;
        for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
          def[i] = warpedPoint[i] - point[i];
        }
        defIter.Set( def );
      }
      this->GetDisplacementField()->Modified();

      typename FFTWRealToComplexType::Pointer fft_filter = FFTWRealToComplexType::New();
      fft_filter->SetInput( this->GetDisplacementField() );
      fft_filter->Update();

      ImageRegionIterator< DisplacementFieldFFTType > coeffsIter( m_Coefficient, m_Coefficient->GetRequestedRegion() );
      ImageRegionConstIterator< ComplexImageType > smoothIter( m_SmoothFilter, m_SmoothFilter->GetRequestedRegion() );
      ImageRegionConstIterator< DisplacementFieldFFTType > iter( fft_filter->GetOutput(),
                                                                 fft_filter->GetOutput()->GetRequestedRegion() );
      for ( iter.GoToBegin(), smoothIter.GoToBegin(), coeffsIter.GoToBegin(); !iter.IsAtEnd();
            ++iter, ++smoothIter, ++coeffsIter )
      {
        typename DisplacementFieldFFTType::PixelType pixel1 = coeffsIter.Get();
        for ( unsigned int i = 0; i < 3; i++ )
        {
          pixel1[i] = pixel1[i] - ( m_LandmarkWeight * smoothIter.Get() * iter.Get()[i] );
        }
        coeffsIter.Set( pixel1 );
      }
      m_Coefficient->Modified();
    }
#  endif
  }

  typename FFTWComplexToRealType::Pointer inverseFFT = FFTWComplexToRealType::New();
  inverseFFT->SetInput( m_Coefficient );
  inverseFFT->Update();
  m_UpdateBuffer = inverseFFT->GetOutput();
  m_UpdateBuffer->DisconnectPipeline();

  // compute the inverse deformation field
  using InverseDisplacementFieldImageType =
    ICCIterativeInverseDisplacementFieldImageFilter< TDisplacementField, TDisplacementField >;
  typename InverseDisplacementFieldImageType::Pointer inverse = InverseDisplacementFieldImageType::New();
  inverse->SetInput( m_UpdateBuffer );
  inverse->SetStopValue( 1.0e-6 );
  inverse->SetNumberOfIterations( 100 );
  inverse->Update();
  m_InverseUpdateBuffer = inverse->GetOutput();
  m_InverseUpdateBuffer->DisconnectPipeline();
}

/**
 * Compute update at a non boundary neighbourhood
 */

template < typename TFixedImage, typename TMovingImage, typename TDisplacementField >
void
ICCDeformableFunction< TFixedImage, TMovingImage, TDisplacementField >::ComputeMetric( void * gd )
{
  GlobalDataStruct * globalData = reinterpret_cast< GlobalDataStruct * >( gd );

  // WARNING!! We compute the global data without taking into account the current update step.
  // There are several reasons for that: If an exponential, a smoothing or any other operation
  // is applied on the update field, we cannot compute the newMappedCenterPoint here; and even
  // if we could, this would be an often unnecessary time-consuming task.

  // using ConstIterationImageType = ImageRegionConstIterator<TFixedImage>;
  ConstIterationImageType it_f( this->GetFixedImage(), this->GetFixedImage()->GetRequestedRegion() );
  // using IterationImageType = ImageRegionIterator<TFixedImage>;
  IterationImageType it_m( this->m_MovingImageWarper->GetOutput(),
                           this->m_MovingImageWarper->GetOutput()->GetRequestedRegion() );

  for ( it_f.GoToBegin(), it_m.GoToBegin(); !it_f.IsAtEnd(); ++it_f, ++it_m )
  {
    if ( this->GetMovingImageMask() || this->GetFixedImageMask() )
    {
      const double fixedValue = static_cast< double >( it_f.Get() );

      MovingPixelType movingPixValue = static_cast< double >( it_m.Get() );

      const double movingValue = static_cast< double >( movingPixValue );

      PixelType update;
      update.Fill( 0.0 );

      const double speedValue = fixedValue - movingValue;

      typename TFixedImage::IndexType    index = it_f.GetIndex();
      typename FixedImageType::PointType fixedPoint;
      this->GetFixedImage()->TransformIndexToPhysicalPoint( index, fixedPoint );
      if ( this->GetFixedImageMask()->IsInsideInWorldSpace( fixedPoint ) ||
           ( m_MovingMaskImageWarper->GetOutput()->GetPixel( index ) > 0 ) )
      {
        if ( globalData )
        {
          globalData->m_SumOfSquaredDifference += itk::Math::sqr( speedValue );
          globalData->m_NumberOfPixelsProcessed += 1;
          //    globalData->m_SumOfSquaredChange += update.GetSquaredNorm();
        }
      }
    }
    else
    {
      const double fixedValue = static_cast< double >( it_f.Get() );

      MovingPixelType movingPixValue = static_cast< double >( it_m.Get() );

      const double movingValue = static_cast< double >( movingPixValue );

      PixelType update;
      update.Fill( 0.0 );

      const double speedValue = fixedValue - movingValue;

      //    std::cout<<"speedValue"<<speedValue<<std::endl;

      if ( globalData )
      {
        globalData->m_SumOfSquaredDifference += itk::Math::sqr( speedValue );
        globalData->m_NumberOfPixelsProcessed += 1;
        //    globalData->m_SumOfSquaredChange += update.GetSquaredNorm();
      }
    }
  }
  //   std::cout<<"SumOfSquaredDifference"<<globalData->m_SumOfSquaredDifference<<std::endl;
}

template < typename TFixedImage, typename TMovingImage, typename TDisplacementField >
typename ICCDeformableFunction< TFixedImage, TMovingImage, TDisplacementField >::PixelType
ICCDeformableFunction< TFixedImage, TMovingImage, TDisplacementField >::ComputeUpdate(
  const NeighborhoodType & /*it*/, void * /*gd */, const FloatOffsetType & itkNotUsed( offset ) )
{
  PixelType update;

  update.Fill( 0.0 );
  return update;
}

/**
 * Update the metric and release the per-thread-global data.
 */
template < typename TFixedImage, typename TMovingImage, typename TDisplacementField >
void
ICCDeformableFunction< TFixedImage, TMovingImage, TDisplacementField >::ReleaseGlobalDataPointer( void * gd ) const
{
  GlobalDataStruct * globalData = reinterpret_cast< GlobalDataStruct * >( gd );

  m_MetricCalculationLock.lock();
  m_SumOfSquaredDifference += globalData->m_SumOfSquaredDifference;
  m_NumberOfPixelsProcessed += globalData->m_NumberOfPixelsProcessed;
  m_SumOfSquaredChange += globalData->m_SumOfSquaredChange;
  if ( m_NumberOfPixelsProcessed )
  {
    m_Metric = m_SumOfSquaredDifference / static_cast< double >( m_NumberOfPixelsProcessed );
    m_RMSChange = std::sqrt( m_SumOfSquaredChange / static_cast< double >( m_NumberOfPixelsProcessed ) );
  }
  m_MetricCalculationLock.unlock();

  delete globalData;
}
} // end namespace itk

#endif
