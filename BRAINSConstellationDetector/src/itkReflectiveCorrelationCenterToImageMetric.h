/*
 * Author: Hans J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

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
#ifndef __itkReflectiveCorrelationCenterToImageMetric_h
#define __itkReflectiveCorrelationCenterToImageMetric_h

#include "itkPoint.h"
#include "itkIdentityTransform.h"

#include "itkPowellOptimizerv4.h"
#include "landmarksConstellationCommon.h"
#include "GenericTransformImage.h"
#include "itkStatisticsImageFilter.h"
#include "itkNumberToString.h"
#include "itkCompensatedSummation.h"

// Optimize the A,B,C vector
template<typename TOptimizerType>
class Rigid3DCenterReflectorFunctor : public itk::ObjectToObjectMetricBase
{
public:
  typedef Rigid3DCenterReflectorFunctor   Self;
  typedef itk::ObjectToObjectMetricBase   Superclass;
  typedef itk::SmartPointer<Self>         Pointer;
  typedef itk::SmartPointer<const Self>   ConstPointer;
  itkNewMacro( Self );

  enum { SpaceDimension=3 };

  typedef Superclass::ParametersType      ParametersType;
  typedef Superclass::DerivativeType      DerivativeType;
  typedef Superclass::MeasureType         MeasureType;

  typedef TOptimizerType                        OptimizerType;
  typedef typename OptimizerType::Pointer       OptimizerPointer;
  typedef itk::CompensatedSummation< double >   CompensatedSummationType;

  Rigid3DCenterReflectorFunctor() :
  m_params(),
  m_OriginalImage(ITK_NULLPTR),
  m_ResamplerReferenceImage(ITK_NULLPTR),
  m_CenterOfHeadMass(),
  m_BackgroundValue(0),
  m_CenterOfImagePoint(),
  m_Translation(),
  m_imInterp(ITK_NULLPTR),
  m_cc(0.0),
  m_HasLocalSupport(false)
  {
    this->m_Optimizer = OptimizerType::New();
    this->m_Optimizer->SetMetric( &( *this ) );
    this->m_Optimizer->SetStepLength( 0.075 );
    this->m_Optimizer->SetStepTolerance( 1e-8 );
    this->m_Optimizer->SetValueTolerance( 1e-8 );
    this->m_Optimizer->SetMaximumIteration( 50 );

    this->m_params.set_size(SpaceDimension);
    this->m_params.fill(0.0);

    this->m_imInterp = LinearInterpolatorType::New();
  }

  ////////////////////////
  // Mandatory metric functions
  virtual void Initialize(void) throw ( itk::ExceptionObject )
  {
    ParametersType params;
    params.set_size(SpaceDimension);

    params[0] = 0;
    params[1] = 0;
    params[2] = 0;
    // Initialize with current guess;
    this->m_params = params;
    double max_cc = this->GetValue();
    const double HARange = 25.0;
    const double BARange = 15.0;

    // rough search in neighborhood.
    const double one_degree = 1.0F * vnl_math::pi / 180.0F;
    const double HAStepSize = HARange * one_degree * .1;
    const double BAStepSize = BARange * one_degree * .1;
    // Let the powell optimizer do all the work for determining the proper
    // offset
    // Quick search just needs to get an approximate angle correct.
      {
      for( double HA = -HARange * one_degree; HA <= HARange * one_degree; HA += HAStepSize )
        {
        for( double BA = -BARange * one_degree; BA <= BARange * one_degree; BA += BAStepSize )
          {
          const double Offset = 0.0;
          params[0] = HA;
          params[1] = BA;
          params[2] = Offset;

          const double current_cc = this->f(params);
          if( current_cc < max_cc )
            {
            this->m_params = params;
            max_cc = current_cc;
            }
          }
        }
      }
    // DEBUGGING INFORMATION
    if( LMC::globalverboseFlag == true )
      {
      itk::NumberToString<double> doubleToString;
      std::cout << "quick search 15 deg "
                << " HA= " << doubleToString(this->m_params[0] * 180.0 / vnl_math::pi)
                << " BA= " << doubleToString(this->m_params[1] * 180.0 / vnl_math::pi)
                << " XO= " << doubleToString(this->m_params[2])
                << " cc="  <<  doubleToString(this->GetValue())
                << " iterations=" << this->m_Optimizer->GetCurrentIteration()
                << std::endl;
      }
  }

  double f(const ParametersType & params) const
  {
    const double        MaxUnpenalizedAllowedDistance = 8.0;
    const double        DistanceFromCenterOfMass = vcl_abs(params[2]);
    static const double FortyFiveDegreesAsRadians = 45.0 * vnl_math::pi / 180.0;
    const double        cost_of_HeadingAngle = ( vcl_abs(params[0]) < FortyFiveDegreesAsRadians ) ? 0 :
      ( ( vcl_abs(params[0]) - FortyFiveDegreesAsRadians ) * 2 );
    const double cost_of_BankAngle = ( vcl_abs(params[1]) < FortyFiveDegreesAsRadians ) ? 0 :
      ( ( vcl_abs(params[1]) - FortyFiveDegreesAsRadians ) * 2 );

    if( ( vcl_abs(params[0]) > FortyFiveDegreesAsRadians ) || ( vcl_abs(params[1]) > FortyFiveDegreesAsRadians ) )
      {
      std::cout << "WARNING: ESTIMATED ROTATIONS ARE WAY TOO BIG SO GIVING A HIGH COST" << std::endl;
      return 1;
      }
    const double cc = -CenterImageReflection_crossCorrelation(params);

    const double cost_of_motion = ( vcl_abs(DistanceFromCenterOfMass) < MaxUnpenalizedAllowedDistance ) ? 0 :
      ( vcl_abs(DistanceFromCenterOfMass - MaxUnpenalizedAllowedDistance) * .1 );
    const double raw_finalcos_gamma = cc + cost_of_motion + cost_of_BankAngle + cost_of_HeadingAngle;

#ifdef __USE_EXTENSIVE_DEBUGGING__
    if( !vnl_math_isfinite(raw_finalcos_gamma) )
      {
      std::cout << __FILE__ << " " << __LINE__ << " "
                << params << " : " << cc << " " << cost_of_HeadingAngle << " " << cost_of_BankAngle << " "
                << cost_of_motion << std::endl;
      return EXIT_FAILURE;
      }
#endif
    return raw_finalcos_gamma;
  }

  virtual MeasureType GetValue() const
  {
    return f(this->m_params);
  }

  virtual void GetDerivative( DerivativeType & ) const
  {
  }

  void GetValueAndDerivative( MeasureType & value,
                              DerivativeType & derivative ) const
  {
    value = GetValue();
    GetDerivative( derivative );
  }

  virtual unsigned int GetNumberOfLocalParameters() const
  {
    return SpaceDimension;
  }

  virtual unsigned int GetNumberOfParameters(void) const
  {
    return SpaceDimension;
  }

  virtual void SetParameters( ParametersType & parameters )
  {
    m_params = parameters;
  }

  virtual const ParametersType & GetParameters() const
  {
    return m_params;
  }

  virtual bool HasLocalSupport() const
  {
    return m_HasLocalSupport;
  }

  void SetHasLocalSupport(bool hls)
  {
    m_HasLocalSupport = hls;
  }

  virtual void UpdateTransformParameters( const DerivativeType &, ParametersValueType )
  {
  }
  ////////////////////////

  void SetCenterOfHeadMass(SImageType::PointType centerOfHeadMass)
  {
    m_CenterOfHeadMass = centerOfHeadMass;
  }

  RigidTransformType::Pointer GetTransformToMSP(void) const
  {
    // Compute and store the new output image origin
    SImageType::Pointer image = GetResampledImageToOutputBox(this->m_params);

    // it is also the msp location
    SImageType::PointType physCenter = GetImageCenterPhysicalPoint(image);

    // Move the physical origin to the center of the image
    RigidTransformType::Pointer tempEulerAngles3DT = RigidTransformType::New();
    tempEulerAngles3DT->Compose( this->GetTransformFromParams(this->m_params) );
    RigidTransformType::TranslationType tnsl = tempEulerAngles3DT->GetTranslation();

    tempEulerAngles3DT->Translate(physCenter.GetVectorFromOrigin() - tnsl);
    return tempEulerAngles3DT;
  }

  void InitializeImage(SImageType::Pointer & RefImage)
  {
      {
      SImageType::PixelType dummy;
      // Find threshold below which image is considered background.
      m_BackgroundValue = setLowHigh<SImageType>(RefImage, dummy, dummy, 0.1F);
      }

    this->m_CenterOfImagePoint = GetImageCenterPhysicalPoint(RefImage);

#ifdef USE_DEBUGGIN_IMAGES
      {
      itk::ImageFileWriter<SImageType>::Pointer dbgWriter = itk::ImageFileWriter<SImageType>::New();
      dbgWriter->UseCompressionOn();
      dbgWriter->SetFileName("itkReflectiveCorrelationCenterToImageMetric149.nii.gz");
      dbgWriter->SetInput(RefImage);
      dbgWriter->Update();
      }
#endif

    this->m_OriginalImage = RefImage;
    // Update the output reference image for the resampler every time the OriginalImage is updated
    this->CreateResamplerReferenceImage();

    this->m_Translation = this->m_CenterOfHeadMass.GetVectorFromOrigin() - m_CenterOfImagePoint.GetVectorFromOrigin();
    if( LMC::globalverboseFlag == true )
      {
      std::cout << "Center Of Physical Point: " << this->m_CenterOfImagePoint << std::endl;
      std::cout << "Center Of Mass Point:" << this->m_CenterOfHeadMass << std::endl;
      std::cout << "InitialTranslation: " << this->m_Translation << std::endl;
      }
  }

  /* -- */
  void SetDownSampledReferenceImage(SImageType::Pointer & NewImage)
  {
    this->m_OriginalImage = NewImage;
    // Update the output reference image for the resampler every time the OriginalImage is updated
    this->CreateResamplerReferenceImage();
  }

  /* -- */
  SImageType::PixelType GetBackgroundValue(void) const
  {
    return this->m_BackgroundValue;
  }

  /* -- */
  RigidTransformType::Pointer GetTransformFromParams(ParametersType const & params) const
  {
    RigidTransformType::Pointer tempEulerAngles3DT = RigidTransformType::New();

    tempEulerAngles3DT->SetCenter(this->m_CenterOfImagePoint);
    tempEulerAngles3DT->SetRotation(0, params[1], params[0]);
    SImageType::PointType::VectorType tnsl = this->m_Translation;
    tnsl[0] += params[2];
    tempEulerAngles3DT->SetTranslation(tnsl);
    return tempEulerAngles3DT;
  }

  void CreateResamplerReferenceImage(void)
  {
    SImageType::SpacingType           outputImageSpacing;
    SImageType::SizeType              outputImageSize;
    SImageType::PointType             outputImageOrigin;
      {
      // Get output spacing
      const SImageType::SpacingType &inputImageSpacing = m_OriginalImage->GetSpacing();
      SImageType::SpacingType::ValueType minSpacing=inputImageSpacing[0];
      for( unsigned int i = 1; i < 3; ++i )
        {
        minSpacing = std::min(minSpacing,inputImageSpacing[i]);
        }
      for( unsigned int i = 0; i < 3; ++i )
        {
        outputImageSpacing[i]=minSpacing;
        }

      // Desire a 95*2 x 130*2 x 160x2 mm voxel lattice that will fit a brain
      outputImageSize[0] = static_cast<unsigned long int>( 2.0 * vcl_ceil(95.0  / outputImageSpacing[0]) );
      outputImageSize[1] = static_cast<unsigned long int>( 2.0 * vcl_ceil(130.0 / outputImageSpacing[1]) );
      outputImageSize[2] = static_cast<unsigned long int>( 2.0 * vcl_ceil(160.0 / outputImageSpacing[2]) );

      // The physical center of MSP plane is not determined yet. At the
      // optimizing stage we take COM as physical center
      outputImageOrigin[0] = m_CenterOfHeadMass[0] - .5 * ( outputImageSize[0] - 1 ) * outputImageSpacing[0];
      outputImageOrigin[1] = m_CenterOfHeadMass[1] - .5 * ( outputImageSize[1] - 1 ) * outputImageSpacing[1];
      outputImageOrigin[2] = m_CenterOfHeadMass[2] - .5 * ( outputImageSize[2] - 1 ) * outputImageSpacing[2];
      }

    // Define start index
    SImageType::IndexType             outputImageStartIndex;
    outputImageStartIndex.Fill(0);

    // Define the output image direction identical
    SImageType::DirectionType         outputImageDirection;
    outputImageDirection.SetIdentity();

    // Define image region
    SImageType::RegionType outputImageRegion;
    outputImageRegion.SetSize(outputImageSize);
    outputImageRegion.SetIndex(outputImageStartIndex);

    this->m_ResamplerReferenceImage = SImageType::New();
    this->m_ResamplerReferenceImage->SetOrigin(outputImageOrigin);
    this->m_ResamplerReferenceImage->SetDirection(outputImageDirection);
    this->m_ResamplerReferenceImage->SetSpacing(outputImageSpacing);
    this->m_ResamplerReferenceImage->SetRegions(outputImageRegion);
    this->m_ResamplerReferenceImage->Allocate();
  }

  /* -- */
  SImageType::Pointer GetResampledImageToOutputBox(ParametersType const & params) const
  {
    /*
     * Resample the image
     */
    ResampleFilterType::Pointer       resampleFilter = ResampleFilterType::New();
    resampleFilter->SetInterpolator(this->m_imInterp);
    resampleFilter->SetDefaultPixelValue(0);
    resampleFilter->UseReferenceImageOn();
    resampleFilter->SetReferenceImage(this->m_ResamplerReferenceImage);
    resampleFilter->SetInput(this->m_OriginalImage);
    resampleFilter->SetTransform( this->GetTransformFromParams(params) );
    resampleFilter->Update();
    return resampleFilter->GetOutput();
  }

  double CenterImageReflection_crossCorrelation(ParametersType const & params) const
  {
    SImageType::Pointer  internalResampledForReflectiveComputationImage = GetResampledImageToOutputBox(params);

    /*
     * Compute the reflective correlation
     */
    double               sumVoxelValues = 0.0F;
    double               sumSquaredVoxelValues = 0.0F;
    double               sumVoxelValuesQR = 0.0F;
    double               sumVoxelValuesReflected = 0.0F;
    double               sumSquaredVoxelValuesReflected = 0.0F;
    int                  N = 0;
    SImageType::SizeType rasterResampleSize = internalResampledForReflectiveComputationImage->GetLargestPossibleRegion().GetSize();
    const SImageType::SizeType::SizeValueType xMaxIndexResampleSize = rasterResampleSize[0] - 1;
    rasterResampleSize[0] /= 2; // Only need to do 1/2 in the x direction;
    SImageType::RegionType rasterRegion;
    rasterRegion.SetSize(rasterResampleSize);
    rasterRegion.SetIndex( internalResampledForReflectiveComputationImage->GetLargestPossibleRegion().GetIndex() );
    itk::ImageRegionConstIteratorWithIndex<SImageType> halfIt(internalResampledForReflectiveComputationImage,
                                                              rasterRegion);

    CompensatedSummationType  CS_sumVoxelValuesQR;
    CompensatedSummationType  CS_sumSquaredVoxelValuesReflected;
    CompensatedSummationType  CS_sumVoxelValuesReflected;
    CompensatedSummationType  CS_sumSquaredVoxelValues;
    CompensatedSummationType  CS_sumVoxelValues;

    for( halfIt.GoToBegin(); !halfIt.IsAtEnd(); ++halfIt )
      {
      // NOTE:  Only need to compute left half of space because of reflection.
      const double _f = halfIt.Get();
      if( _f < this->m_BackgroundValue )  // don't worry about background
                                          // voxels.
        {
        continue;
        }
      SImageType::IndexType ReflectedIndex = halfIt.GetIndex();
      ReflectedIndex[0] = xMaxIndexResampleSize - ReflectedIndex[0];
      const double g = internalResampledForReflectiveComputationImage->GetPixel(ReflectedIndex);
      if( g < this->m_BackgroundValue )  // don't worry about background voxels.
        {
        continue;
        }
      CS_sumVoxelValuesQR += _f * g;
      CS_sumSquaredVoxelValuesReflected += g * g;
      CS_sumVoxelValuesReflected += g;
      CS_sumSquaredVoxelValues += _f * _f;
      CS_sumVoxelValues += _f;
      N++;
      }
    sumVoxelValuesQR = CS_sumVoxelValuesQR.GetSum();
    sumSquaredVoxelValuesReflected = CS_sumSquaredVoxelValuesReflected.GetSum();
    sumVoxelValuesReflected = CS_sumVoxelValuesReflected.GetSum();
    sumSquaredVoxelValues = CS_sumSquaredVoxelValues.GetSum();
    sumVoxelValues = CS_sumVoxelValues.GetSum();

    // ///////////////////////////////////////////////
    if( N == 0
        || ( ( sumSquaredVoxelValues - sumVoxelValues * sumVoxelValues
               / N ) * ( sumSquaredVoxelValuesReflected - sumVoxelValuesReflected * sumVoxelValuesReflected / N ) ) ==
        0.0 )
      {
      return 0.0;
      }
    const double cc =
      ( ( sumVoxelValuesQR - sumVoxelValuesReflected * sumVoxelValues
          / N )
        / vcl_sqrt( ( sumSquaredVoxelValues - sumVoxelValues * sumVoxelValues
                      / N )
                    * ( sumSquaredVoxelValuesReflected - sumVoxelValuesReflected * sumVoxelValuesReflected / N ) ) );
    return cc;
  }

  SImageType::Pointer GetMSPCenteredImage(void)
  {
    // Note GetTransformToMSP() aligns physical origin to the center of MSP.
    // SimpleResampleImage() further aligns physical center of image to physical
    // origin.

    // RigidTransformType::Pointer tempEulerAngles3DT = GetTransformToMSP();
    // return SimpleResampleImage( this->m_OriginalImage, tempEulerAngles3DT );

    // create a colormap lookup table
    typedef itk::StatisticsImageFilter<SImageType> StatisticsFilterType;
    StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
    statisticsFilter->SetInput(this->m_OriginalImage);
    statisticsFilter->Update();
    SImageType::PixelType minPixelValue = statisticsFilter->GetMinimum();

    return TransformResample<SImageType, SImageType>
             ( this->m_OriginalImage.GetPointer(),
             MakeIsoTropicReferenceImage().GetPointer(),
             minPixelValue,
             GetInterpolatorFromString<SImageType>("Linear").GetPointer(),
             GetTransformToMSP().GetPointer() );
  }

  double GetCC(void) const
  {
    return this->m_cc;
  }

  void Update(void)
  {
    itk::NumberToString<double> doubleToString;

    try
      {
      this->m_Optimizer->StartOptimization();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cout << "Exception thrown ! " << std::endl;
      std::cout << "An error occurred during Optimization" << std::endl;
      std::cout << "Location    = " << e.GetLocation()    << std::endl;
      std::cout << "Description = " << e.GetDescription() << std::endl;
      //return EXIT_FAILURE;
      }

    this->m_params = this->m_Optimizer->GetCurrentPosition();
    this->m_cc = this->GetValue();

    std::cout << doubleToString(this->m_params[0] * 180.0 / vnl_math::pi) << " "
              << doubleToString(this->m_params[1] * 180.0 / vnl_math::pi) << " "
              << this->m_params[2] << " cc= "
              << doubleToString(m_cc) << " iters= " << this->m_Optimizer->GetCurrentIteration()
              << std::endl;
  }

private:
  typedef itk::ResampleImageFilter<SImageType, SImageType> ResampleFilterType;

  ParametersType                    m_params;
  SImageType::Pointer               m_OriginalImage;
  SImageType::Pointer               m_ResamplerReferenceImage;
  SImageType::PointType             m_CenterOfHeadMass;
  SImageType::PixelType             m_BackgroundValue;
  SImageType::PointType             m_CenterOfImagePoint;
  SImageType::PointType::VectorType m_Translation;
  OptimizerPointer                  m_Optimizer;
  LinearInterpolatorType::Pointer   m_imInterp;
  double                            m_cc;
  bool                              m_HasLocalSupport;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkReflectiveCorrelationCenterToImageMetric.hxx"
#endif

#endif
