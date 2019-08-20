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
#ifndef _GenericTransformImage_hxx_
#define _GenericTransformImage_hxx_

#include <iostream>
#include "GenericTransformImage.h"
#include "itkResampleInPlaceImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkIO.h"

template < typename InputImageType, typename OutputImageType >
typename OutputImageType::Pointer
TransformResample(
  typename InputImageType::ConstPointer                                   inputImage,
  typename itk::ImageBase< InputImageType::ImageDimension >::ConstPointer ReferenceImage,
  const typename InputImageType::PixelType                                defaultValue,
  typename itk::InterpolateImageFunction<
    InputImageType, typename itk::NumericTraits< typename InputImageType::PixelType >::RealType >::Pointer interp,
  typename itk::Transform< double, 3, 3 >::ConstPointer                                                    transform )
{
  using ResampleImageFilter = typename itk::ResampleImageFilter< InputImageType, OutputImageType >;
  typename ResampleImageFilter::Pointer resample = ResampleImageFilter::New();
  resample->SetInput( inputImage );
  resample->SetTransform( transform.GetPointer() );
  resample->SetInterpolator( interp.GetPointer() );

  if ( ReferenceImage.IsNotNull() )
  {
    resample->SetOutputParametersFromImage( ReferenceImage );
  }
  else
  {
    std::cout << "Alert:  missing Reference Volume information default image size set to inputImage" << std::endl;
    resample->SetOutputParametersFromImage( inputImage );
  }
  resample->SetDefaultPixelValue( defaultValue );
  resample->Update();
  typename OutputImageType::Pointer returnval = resample->GetOutput();

  // returnval->DisconnectPipeline();
  return returnval;
}

template < typename InputImageType, typename OutputImageType, typename DisplacementImageType >
typename OutputImageType::Pointer
TransformWarp(
  InputImageType const * const inputImage, const itk::ImageBase< InputImageType::ImageDimension > * ReferenceImage,
  typename InputImageType::PixelType defaultValue,
  typename itk::InterpolateImageFunction<
    InputImageType, typename itk::NumericTraits< typename InputImageType::PixelType >::RealType >::Pointer interp,
  typename DisplacementImageType::Pointer displacementField )
{
  using WarpImageFilter = typename itk::WarpImageFilter< InputImageType, OutputImageType, DisplacementImageType >;
  typename WarpImageFilter::Pointer warp = WarpImageFilter::New();
  warp->SetInput( inputImage );
  warp->SetDisplacementField( displacementField );
  warp->SetInterpolator( interp );

  if ( ReferenceImage != nullptr )
  {
    warp->SetOutputParametersFromImage( ReferenceImage );
  }
  else
  {
    std::cout << "Alert:  missing Reference Volume information default image size set to displacementField"
              << std::endl;
    warp->SetOutputParametersFromImage( displacementField );
  }
  warp->SetEdgePaddingValue( defaultValue );
  warp->Update();
  typename OutputImageType::Pointer returnval = warp->GetOutput();
  // returnval->DisconnectPipeline();
  return returnval;
}

template < typename InputImageType >
typename itk::InterpolateImageFunction<
  InputImageType, typename itk::NumericTraits< typename InputImageType::PixelType >::RealType >::Pointer
GetInterpolatorFromString( const std::string & interpolationMode )
{
  using TInterpolatorPrecisionType = typename itk::NumericTraits< typename InputImageType::PixelType >::RealType;
  using BoundaryConditionType = typename itk::ConstantBoundaryCondition< InputImageType >;
  //
  // WindowRadius set to 3 to match default for Slicer
  static constexpr unsigned int SlicerWindowedSincWindowRadius = 3;

  if ( interpolationMode == "NearestNeighbor" )
  {
    using InterpolatorType =
      typename itk::NearestNeighborInterpolateImageFunction< InputImageType, TInterpolatorPrecisionType >;
    return ( InterpolatorType::New() ).GetPointer();
  }
  else if ( interpolationMode == "Linear" )
  {
    using InterpolatorType = typename itk::LinearInterpolateImageFunction< InputImageType, TInterpolatorPrecisionType >;
    return ( InterpolatorType::New() ).GetPointer();
  }
  else if ( interpolationMode == "BSpline" )
  {
    using InterpolatorType =
      typename itk::BSplineInterpolateImageFunction< InputImageType, TInterpolatorPrecisionType >;
    return ( InterpolatorType::New() ).GetPointer();
  }
  else if ( interpolationMode == "WindowedSinc" )
  {
    static constexpr unsigned int WindowedSincHammingWindowRadius = 5;
    using WindowFunctionType = typename itk::Function::
      HammingWindowFunction< WindowedSincHammingWindowRadius, TInterpolatorPrecisionType, TInterpolatorPrecisionType >;
    using InterpolatorType = typename itk::WindowedSincInterpolateImageFunction< InputImageType,
                                                                                 WindowedSincHammingWindowRadius,
                                                                                 WindowFunctionType,
                                                                                 BoundaryConditionType,
                                                                                 TInterpolatorPrecisionType >;
    return ( InterpolatorType::New() ).GetPointer();
  }
  else if ( interpolationMode == "Hamming" )
  {
    using WindowFunctionType = typename itk::Function::
      HammingWindowFunction< SlicerWindowedSincWindowRadius, TInterpolatorPrecisionType, TInterpolatorPrecisionType >;
    using InterpolatorType = typename itk::WindowedSincInterpolateImageFunction< InputImageType,
                                                                                 SlicerWindowedSincWindowRadius,
                                                                                 WindowFunctionType,
                                                                                 BoundaryConditionType,
                                                                                 TInterpolatorPrecisionType >;
    return ( InterpolatorType::New() ).GetPointer();
  }
  else if ( interpolationMode == "Cosine" )
  {
    using WindowFunctionType = typename itk::Function::
      CosineWindowFunction< SlicerWindowedSincWindowRadius, TInterpolatorPrecisionType, TInterpolatorPrecisionType >;
    using InterpolatorType = typename itk::WindowedSincInterpolateImageFunction< InputImageType,
                                                                                 SlicerWindowedSincWindowRadius,
                                                                                 WindowFunctionType,
                                                                                 BoundaryConditionType,
                                                                                 TInterpolatorPrecisionType >;
    return ( InterpolatorType::New() ).GetPointer();
  }
  else if ( interpolationMode == "Welch" )
  {
    using WindowFunctionType = typename itk::Function::
      WelchWindowFunction< SlicerWindowedSincWindowRadius, TInterpolatorPrecisionType, TInterpolatorPrecisionType >;
    using InterpolatorType = typename itk::WindowedSincInterpolateImageFunction< InputImageType,
                                                                                 SlicerWindowedSincWindowRadius,
                                                                                 WindowFunctionType,
                                                                                 BoundaryConditionType,
                                                                                 TInterpolatorPrecisionType >;
    return ( InterpolatorType::New() ).GetPointer();
  }
  else if ( interpolationMode == "Lanczos" )
  {
    using WindowFunctionType = typename itk::Function::
      LanczosWindowFunction< SlicerWindowedSincWindowRadius, TInterpolatorPrecisionType, TInterpolatorPrecisionType >;
    using InterpolatorType = typename itk::WindowedSincInterpolateImageFunction< InputImageType,
                                                                                 SlicerWindowedSincWindowRadius,
                                                                                 WindowFunctionType,
                                                                                 BoundaryConditionType,
                                                                                 TInterpolatorPrecisionType >;
    return ( InterpolatorType::New() ).GetPointer();
  }
  else if ( interpolationMode == "Blackman" )
  {
    using WindowFunctionType = typename itk::Function::
      BlackmanWindowFunction< SlicerWindowedSincWindowRadius, TInterpolatorPrecisionType, TInterpolatorPrecisionType >;
    using InterpolatorType = typename itk::WindowedSincInterpolateImageFunction< InputImageType,
                                                                                 SlicerWindowedSincWindowRadius,
                                                                                 WindowFunctionType,
                                                                                 BoundaryConditionType,
                                                                                 TInterpolatorPrecisionType >;
    return ( InterpolatorType::New() ).GetPointer();
  }
  else
  {
    std::cout << "Error: Invalid interpolation mode specified -" << interpolationMode << "- " << std::endl;
    std::cout << "\tValid modes: NearestNeighbor, Linear, BSpline, WindowedSinc" << std::endl;
  }
  return nullptr;
}

template < typename InputImageType, typename OutputImageType, typename DisplacementImageType >
typename OutputImageType::Pointer
GenericTransformImage( InputImageType const * const                             OperandImage,
                       const itk::ImageBase< InputImageType::ImageDimension > * ReferenceImage,
                       // typename DisplacementImageType::Pointer DisplacementField,
                       typename itk::Transform< double, 3, 3 >::ConstPointer genericTransform,
                       typename InputImageType::PixelType                    suggestedDefaultValue, // NOTE:  This is
                                                                                                    // ignored in the
                                                                                                    // case of binary
                                                                                                    // image!
                       const std::string & interpolationMode, const bool binaryFlag )
{
  // FIRST will need to convert binary image to signed distance in case
  // binaryFlag is true.

  typename InputImageType::ConstPointer PrincipalOperandImage;

  // Splice in a case for dealing with binary images,
  // where signed distance maps are warped and thresholds created.
  if ( binaryFlag )
  {
    if ( interpolationMode == "NearestNeighbor" )
    {
      std::cout << "WARNING:  Using NearestNeighbor and SignedDistance" << std::endl
                << "          for binary images is an unlikely combination." << std::endl
                << "          you probably want Linear interpolationMode for" << std::endl
                << "          the signed distance map implied by your choice" << std::endl
                << "          of pixelType binary." << std::endl;
    }
    /* We make the values inside the structures positive and outside negative
     * using
     *  BinaryThresholdImageFilter. As the lower and upper threshold values are
     *     0 only values of 0 in the image are filled with 0.0 and other
     *     values are  1.0
     */

    using FloatThresholdFilterType = itk::BinaryThresholdImageFilter< InputImageType, InputImageType >;
    typename FloatThresholdFilterType::Pointer initialFilter = FloatThresholdFilterType::New();
    initialFilter->SetInput( OperandImage );
    {
      constexpr typename FloatThresholdFilterType::OutputPixelType outsideValue = 1.0;
      constexpr typename FloatThresholdFilterType::OutputPixelType insideValue = 0.0;
      initialFilter->SetOutsideValue( outsideValue );
      initialFilter->SetInsideValue( insideValue );
      constexpr typename FloatThresholdFilterType::InputPixelType lowerThreshold = 0;
      constexpr typename FloatThresholdFilterType::InputPixelType upperThreshold = 0;
      initialFilter->SetLowerThreshold( lowerThreshold );
      initialFilter->SetUpperThreshold( upperThreshold );
    }
    initialFilter->Update();
    {
      using DistanceFilterType = itk::SignedMaurerDistanceMapImageFilter< InputImageType, InputImageType >;
      typename DistanceFilterType::Pointer DistanceFilter = DistanceFilterType::New();
      DistanceFilter->SetInput( initialFilter->GetOutput() );
      // DistanceFilter->SetNarrowBandwidth( m_BandWidth );
      DistanceFilter->SetInsideIsPositive( true );
      DistanceFilter->SetUseImageSpacing( true );
      DistanceFilter->SetSquaredDistance( false );

      DistanceFilter->Update();
      PrincipalOperandImage = DistanceFilter->GetOutput();
      // PrincipalOperandImage->DisconnectPipeline();
    }
    // Using suggestedDefaultValue based on the size of the image so that
    // intensity values
    // are kept to a reasonable range.  (A costlier way calculates the image
    // min.)
    const typename InputImageType::SizeType    size = PrincipalOperandImage->GetLargestPossibleRegion().GetSize();
    const typename InputImageType::SpacingType spacing = PrincipalOperandImage->GetSpacing();
    double                                     diagonalLength = 0;
    for ( unsigned int s = 0; s < InputImageType::ImageDimension; ++s )
    {
      diagonalLength += size[s] * spacing[s];
    }
    // Consider the 3D diagonal value, to guarantee that the background
    // filler is unlikely to add shapes to the thresholded signed
    // distance image. This is an easy enough proof of a lower bound on
    // the image min, since it works even if the mask is a single voxel in
    // the image field corner. suggestedDefaultValue=
    // std::sqrt( diagonalLength );
    // In most cases, a heuristic fraction of the diagonal value is an
    // even better lower bound: if the midpoint of the image is inside the
    // mask, 1/2 is a lower bound as well, and the background is unlikely
    // to drive the upper limit of the intensity range when we visualize
    // the intermediate image for debugging.

    suggestedDefaultValue = -std::sqrt( diagonalLength ) * 0.5;
  }
  else // other than if (pixelType == "binary")
  {
    PrincipalOperandImage = OperandImage;
  }

  // RESAMPLE with the appropriate transform and interpolator:
  // One name for the intermediate resampled float image.
  typename InputImageType::Pointer TransformedImage;

  // std::cout<< " Displacement Field is Null... " << std::endl;
  if ( genericTransform.IsNotNull() )
  {
    if ( interpolationMode == "ResampleInPlace" )
    {
      using ResampleIPFilterType = itk::ResampleInPlaceImageFilter< InputImageType, OutputImageType >;
      using ResampleIPFilterPointer = typename ResampleIPFilterType::Pointer;

      using CompositeTransformType = itk::CompositeTransform< double, 3 >;
      using VersorRigid3DTransformType = itk::VersorRigid3DTransform< double >;

      VersorRigid3DTransformType::Pointer tempInitializerITKTransform = VersorRigid3DTransformType::New();

      std::string genericTransformFileType;
      if ( genericTransform.IsNotNull() )
      {
        genericTransformFileType = genericTransform->GetNameOfClass();
      }
      if ( genericTransformFileType == "CompositeTransform" )
      {
        const CompositeTransformType::ConstPointer genericCompositeTransform =
          static_cast< const CompositeTransformType * >( genericTransform.GetPointer() );
        if ( genericCompositeTransform->GetNumberOfTransforms() > 1 )
        {
          std::cout << "Error in type conversion. " << __FILE__ << __LINE__ << std::endl;
          std::cout << "ResampleInPlace is only allowed with rigid transform type,"
                    << "but the input composite transform consists of more than one transfrom." << std::endl;
        }
        // extract the included linear rigid transform from the input composite
        CompositeTransformType::TransformType::Pointer genericTempTransform =
          genericCompositeTransform->GetNthTransform( 0 ).GetPointer();
        std::string genericTempTransformFileType;
        if ( genericTempTransform.IsNotNull() )
        {
          genericTempTransformFileType = genericTempTransform->GetNameOfClass();
        }
        if ( genericTempTransformFileType == "ScaleSkewVersor3DTransform" ||
             genericTempTransformFileType == "ScaleVersor3DTransform" ||
             genericTempTransformFileType == "Similarity3DTransform" ||
             genericTempTransformFileType == "VersorRigid3DTransform" )
        {
          const VersorRigid3DTransformType::ConstPointer tempTransform =
            static_cast< const VersorRigid3DTransformType * >( genericTempTransform.GetPointer() );
          tempInitializerITKTransform->SetFixedParameters( tempTransform->GetFixedParameters() );
          tempInitializerITKTransform->SetParameters( tempTransform->GetParameters() );
        }
        else
        {
          std::cout << "Error in type conversion. " << __FILE__ << __LINE__ << std::endl;
          std::cout << "ResampleInPlace is only allowed with rigid transform type." << std::endl;
        }
      }
      else
      {
        std::string genericTempTransformFileType;
        if ( genericTransform.IsNotNull() )
        {
          genericTempTransformFileType = genericTransform->GetNameOfClass();
        }
        if ( genericTempTransformFileType == "ScaleSkewVersor3DTransform" ||
             genericTempTransformFileType == "ScaleVersor3DTransform" ||
             genericTempTransformFileType == "Similarity3DTransform" ||
             genericTempTransformFileType == "VersorRigid3DTransform" )
        {
          const VersorRigid3DTransformType::ConstPointer tempTransform =
            static_cast< VersorRigid3DTransformType const * >( genericTransform.GetPointer() );
          tempInitializerITKTransform->SetFixedParameters( tempTransform->GetFixedParameters() );
          tempInitializerITKTransform->SetParameters( tempTransform->GetParameters() );
        }
        else
        {
          std::cout << "Error in type conversion. " << __FILE__ << __LINE__ << std::endl;
          std::cout << "ResampleInPlace is only allowed with rigid transform type." << std::endl;
        }
      }
      if ( tempInitializerITKTransform.IsNull() )
      {
        itkGenericExceptionMacro(
          << "Error in type conversion. ResampleInPlace is only allowed with rigid transform type." );
      }

      ResampleIPFilterPointer resampleIPFilter = ResampleIPFilterType::New();
      resampleIPFilter->SetInputImage( PrincipalOperandImage.GetPointer() );
      resampleIPFilter->SetRigidTransform( tempInitializerITKTransform.GetPointer() );
      resampleIPFilter->Update();
      TransformedImage = resampleIPFilter->GetOutput();
    }
    else
    {
      TransformedImage = TransformResample< InputImageType, OutputImageType >(
        PrincipalOperandImage.GetPointer(),
        // INFO:  Change function signature to be a ConstPointer instead of a raw pointer ReferenceImage.GetPointer(),
        ReferenceImage,
        suggestedDefaultValue,
        GetInterpolatorFromString< InputImageType >( interpolationMode ).GetPointer(),
        genericTransform.GetPointer() );
    }
  }

  // FINALLY will need to convert signed distance to binary image in case
  // binaryFlag is true.

  typename InputImageType::Pointer FinalTransformedImage;

  if ( binaryFlag )
  {
    // A special case for dealing with binary images
    // where signed distance maps are warped and thresholds created
    using MaskPixelType = short int;
    using BinFlagOnMaskImageType = typename itk::Image< MaskPixelType, 3 >;

    // Now Threshold and write out image
    using BinaryThresholdFilterType =
      typename itk::BinaryThresholdImageFilter< InputImageType, BinFlagOnMaskImageType >;
    typename BinaryThresholdFilterType::Pointer finalFilter = BinaryThresholdFilterType::New();
    finalFilter->SetInput( TransformedImage );

    constexpr typename BinaryThresholdFilterType::OutputPixelType outsideValue = 0;
    constexpr typename BinaryThresholdFilterType::OutputPixelType insideValue = 1;
    finalFilter->SetOutsideValue( outsideValue );
    finalFilter->SetInsideValue( insideValue );
    // Signed distance boundary voxels are defined as being included in the
    // structure,  therefore the desired distance threshold is in the middle
    // of the enclosing (negative) voxel ribbon around threshold 0.
    const typename InputImageType::SpacingType               Spacing = ReferenceImage->GetSpacing();
    const typename BinaryThresholdFilterType::InputPixelType lowerThreshold =
      -0.5 * 0.333333333333 * ( Spacing[0] + Spacing[1] + Spacing[2] );
    //  std::cerr << "Lower Threshold == " << lowerThreshold << std::endl;

    const typename BinaryThresholdFilterType::InputPixelType upperThreshold =
      std::numeric_limits< typename BinaryThresholdFilterType::InputPixelType >::max();
    finalFilter->SetLowerThreshold( lowerThreshold );
    finalFilter->SetUpperThreshold( upperThreshold );

    finalFilter->Update();

    using CastImageFilter = typename itk::CastImageFilter< BinFlagOnMaskImageType, InputImageType >;
    typename CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( finalFilter->GetOutput() );
    castFilter->Update();

    FinalTransformedImage = castFilter->GetOutput();
  }
  else
  {
    FinalTransformedImage = TransformedImage;
  }

  return FinalTransformedImage;
}

#endif
