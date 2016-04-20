/*
 * Author: Hans J. Johnson, Wei Lu, Ali Ghayoor
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
  m_CenterOfHeadMassIsSet(false),
  m_BackgroundValue(0),
  m_DoPowell(true),
  m_imInterp(ITK_NULLPTR),
  m_cc(0.0),
  m_HasLocalSupport(false)
  {
    this->m_Optimizer = OptimizerType::New();
    this->m_Optimizer->SetMetric( &( *this ) );
    this->m_Optimizer->SetStepLength( 0.025 );
    this->m_Optimizer->SetStepTolerance( 1e-12 );
    this->m_Optimizer->SetValueTolerance( 1e-12 );
    this->m_Optimizer->SetMaximumIteration( 50 );

    this->m_params.set_size(SpaceDimension);
    this->m_params.fill(0.0);

    this->m_imInterp = LinearInterpolatorType::New();
  }

  ////////////////////////
  // Mandatory metric functions
  virtual void Initialize(void) throw ( itk::ExceptionObject ) ITK_OVERRIDE
  {
    ParametersType params;
    params.set_size(SpaceDimension);

    params[0] = 0;
    params[1] = 0;
    params[2] = 0;
    // Initialize with current guess;
    this->m_params = params;
    double max_cc = this->GetValue();

    // Run a multi-level exhaustive search to find the approximate parameters.
    // The Powell optimizer will then tune the exhaustive search parameters.
    std::vector<double> Angle_Range(3);
    Angle_Range[0] = 45.0;
    Angle_Range[1] = 2.5;
    Angle_Range[2] = 0.5;

    std::vector<double> Angle_Stepsizes(3);
    Angle_Stepsizes[0] = 5.0;
    Angle_Stepsizes[1] = 0.5;
    Angle_Stepsizes[2] = 0.25;

    std::vector<double> Offset_Range(3);
    Offset_Range[0] = 15.0;
    Offset_Range[1] = 1.5;
    Offset_Range[2] = 0.5;

    std::vector<double> Offset_Stepsizes(3);
    Offset_Stepsizes[0] = 3.0;
    Offset_Stepsizes[1] = 0.5;
    Offset_Stepsizes[2] = 0.25;

    for( unsigned int resolutionIter = 0; resolutionIter <= 2; ++resolutionIter )
      {
      const double HA_range = Angle_Range[resolutionIter];
      const double BA_range = Angle_Range[resolutionIter];
      const double LR_range = Offset_Range[resolutionIter];

      const double HA_stepsize = Angle_Stepsizes[resolutionIter]; // degree
      const double BA_stepsize = Angle_Stepsizes[resolutionIter]; // degree
      const double LR_stepsize = Offset_Stepsizes[resolutionIter]; // mm

      this->DoExhaustiveSearch(this->m_params, max_cc,
                               HA_range, BA_range, LR_range,
                               HA_stepsize, BA_stepsize, LR_stepsize
#ifdef WRITE_CSV_FILE
                               ,std::string("")
#endif
                               );
      }

    this->m_cc = max_cc;

    // DEBUGGING INFORMATION
    if( LMC::globalverboseFlag )
      {
      itk::NumberToString<double> doubleToString;
      std::cout << "Initialize exhaustive search: "
                << " HA= " << doubleToString(this->m_params[0] * 180.0 / vnl_math::pi)
                << " BA= " << doubleToString(this->m_params[1] * 180.0 / vnl_math::pi)
                << " XO= " << doubleToString(this->m_params[2])
                << " cc="  <<  doubleToString(this->m_cc)
                << std::endl;
      }
  }

  virtual MeasureType GetValue() const ITK_OVERRIDE
  {
    return f(this->m_params);
  }

  virtual void GetDerivative( DerivativeType & ) const ITK_OVERRIDE
  {
  }

  void GetValueAndDerivative( MeasureType & value,
                              DerivativeType & derivative ) const ITK_OVERRIDE
  {
    value = GetValue();
    GetDerivative( derivative );
  }

  virtual unsigned int GetNumberOfLocalParameters() const ITK_OVERRIDE
  {
    return SpaceDimension;
  }

  virtual unsigned int GetNumberOfParameters(void) const ITK_OVERRIDE
  {
    return SpaceDimension;
  }

  virtual void SetParameters( ParametersType & parameters ) ITK_OVERRIDE
  {
    this->m_params = parameters;
  }

  virtual const ParametersType & GetParameters() const ITK_OVERRIDE
  {
    return this->m_params;
  }

  virtual bool HasLocalSupport() const ITK_OVERRIDE
  {
    return this->m_HasLocalSupport;
  }

  void SetHasLocalSupport(bool hls)
  {
    this->m_HasLocalSupport = hls;
  }

  virtual void UpdateTransformParameters( const DerivativeType &, ParametersValueType ) ITK_OVERRIDE
  {
  }
  ////////////////////////

  void DoExhaustiveSearch(ParametersType & opt_params,
                          double & opt_cc,
                          const double HARange,
                          const double BARange,
                          const double LRRange,
                          const double HAStepSize,
                          const double BAStepSize,
                          const double LRStepSize
#ifdef WRITE_CSV_FILE
                          ,std::string CSVFileName
#endif
                          )
  {
  // starting parameters are optimal parameters from previous level
  const ParametersType starting_params (opt_params);
  // new search parameters
  ParametersType current_params;
  current_params.set_size(SpaceDimension);
#ifdef WRITE_CSV_FILE
  std::stringstream csvFileOfMetricValues;
#endif
  const double degree_to_rad = vnl_math::pi / 180.0;

  for( double LR = -LRRange; LR <= LRRange; LR += LRStepSize)
    {
    for( double HA = -HARange; HA <= HARange; HA += HAStepSize )
      {
      for( double BA = -BARange; BA <= BARange; BA += BAStepSize )
        {
        current_params[0] = starting_params[0]+HA * degree_to_rad;
        current_params[1] = starting_params[1]+BA * degree_to_rad;
        current_params[2] = starting_params[2]+LR;

        const double current_cc = this->f(current_params);
        if( current_cc < opt_cc )
          {
          opt_params = current_params;
          opt_cc = current_cc;
          }

#ifdef WRITE_CSV_FILE
        csvFileOfMetricValues << current_params[0]/degree_to_rad
                              << "," << current_params[1]/degree_to_rad
                              << "," << current_params[2]
                              << "," << current_cc
                              << std::endl;
#endif
        }
      }
    }
#ifdef WRITE_CSV_FILE
  if( CSVFileName != "" )
    {
    std::cout << "\nWriting out metric values in a csv file..." << std::endl;
    std::ofstream csvFile;
    csvFile.open( CSVFileName.c_str() );
    if( !csvFile.is_open() )
      {
      itkGenericExceptionMacro( << "Error: Can't write oputput csv file: " << CSVFileName << "!" << std::endl );
      }
    csvFile << csvFileOfMetricValues.str();
    csvFile.close();
    }
#endif
  }

  double f(const ParametersType & params) const
  {
  const double        MaxUnpenalizedAllowedDistance = 8.0;
  const double        DistanceFromCenterOfMass = std::abs(params[2]);
  static const double FortyFiveDegreesAsRadians = 45.0 * vnl_math::pi / 180.0;
  const double        cost_of_HeadingAngle = ( std::abs(params[0]) < FortyFiveDegreesAsRadians ) ? 0 :
  ( ( std::abs(params[0]) - FortyFiveDegreesAsRadians ) * 2 );
  const double cost_of_BankAngle = ( std::abs(params[1]) < FortyFiveDegreesAsRadians ) ? 0 :
  ( ( std::abs(params[1]) - FortyFiveDegreesAsRadians ) * 2 );

  if( ( std::abs(params[0]) > FortyFiveDegreesAsRadians ) || ( std::abs(params[1]) > FortyFiveDegreesAsRadians ) )
    {
    std::cout << "WARNING: ESTIMATED ROTATIONS ARE WAY TOO BIG SO GIVING A HIGH COST" << std::endl;
    return 1;
    }
  const double cc = -CenterImageReflection_crossCorrelation(params);

  const double cost_of_motion = ( std::abs(DistanceFromCenterOfMass) < MaxUnpenalizedAllowedDistance ) ? 0 :
  ( std::abs(DistanceFromCenterOfMass - MaxUnpenalizedAllowedDistance) * .1 );
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


  RigidTransformType::Pointer GetTransformToMSP(void) const
  {
    // Here we try to make MSP plane as the mid slice of the output image voxel lattice
    SImageType::Pointer image = GetResampledImageToOutputBox(this->m_params);
    // it should be the msp location
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
    if( !m_CenterOfHeadMassIsSet )
      {
      itkGenericExceptionMacro(<< "ERROR: m_CenterOfHeadMass is not set!" << std::endl);
      }
    else
      {
      if( LMC::globalverboseFlag )
        {
        std::cout << "Center of Head Mass: [" << this->m_CenterOfHeadMass[0] << "," << this->m_CenterOfHeadMass[1] << ","
          << this->m_CenterOfHeadMass[2] << "]" << std::endl;
        }
      }

      {
      SImageType::PixelType dummyLow;
      SImageType::PixelType dummyHigh;
      // Find threshold below which image is considered background.
      m_BackgroundValue = setLowHigh<SImageType>(RefImage, dummyLow, dummyHigh, 0.1F);
      }

#ifdef USE_DEBUGGIN_IMAGES
      {
      itk::ImageFileWriter<SImageType>::Pointer dbgWriter = itk::ImageFileWriter<SImageType>::New();
      dbgWriter->UseCompressionOn();
      dbgWriter->SetFileName("itkReflectiveCorrelationCenterToImageMetric149.nii.gz");
      dbgWriter->SetInput(RefImage);
      dbgWriter->Update();
      }
#endif

    this->SetDownSampledReferenceImage(RefImage);
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

    tempEulerAngles3DT->SetCenter(this->GetCenterOfHeadMass());
    tempEulerAngles3DT->SetRotation(0, params[1], params[0]);
    SImageType::PointType::VectorType tnsl;
    tnsl[0] = params[2];
    tnsl[1] = 0;
    tnsl[2] = 0;
    tempEulerAngles3DT->SetTranslation(tnsl);

    return tempEulerAngles3DT;
  }

  void CreateResamplerReferenceImage(void)
  {
    SImageType::SizeType              outputImageSize;
    SImageType::PointType             outputImageOrigin;

    SImageType::SpacingType outputImageSpacing = m_OriginalImage->GetSpacing();
      {
      // Get output spacing
      SImageType::SpacingType::ValueType minSpacing=outputImageSpacing[0];
      for( unsigned int i = 1; i < 3; ++i )
        {
        minSpacing = std::min(minSpacing,outputImageSpacing[i]);
        }
      for( unsigned int i = 0; i < 3; ++i )
        {
        outputImageSpacing[i]=minSpacing;
        }

      // Desire a 95*2 x 130*2 x 160x2 mm voxel lattice that will fit a brain
      outputImageSize[0] = static_cast<SImageType::SpacingType::ValueType>( 2.0 * std::ceil(95.0  / outputImageSpacing[0]) );
      outputImageSize[1] = static_cast<SImageType::SpacingType::ValueType>( 2.0 * std::ceil(130.0 / outputImageSpacing[1]) );
      outputImageSize[2] = static_cast<SImageType::SpacingType::ValueType>( 2.0 * std::ceil(160.0 / outputImageSpacing[2]) );

      // The physical center of MSP plane is not determined yet. At the
      // optimizing stage we take COM as physical center
      outputImageOrigin[0] = this->GetCenterOfHeadMass()[0] - .5 * ( outputImageSize[0] - 1 ) * outputImageSpacing[0];
      outputImageOrigin[1] = this->GetCenterOfHeadMass()[1] - .5 * ( outputImageSize[1] - 1 ) * outputImageSpacing[1];
      outputImageOrigin[2] = this->GetCenterOfHeadMass()[2] - .5 * ( outputImageSize[2] - 1 ) * outputImageSpacing[2];
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
        / std::sqrt( ( sumSquaredVoxelValues - sumVoxelValues * sumVoxelValues
                      / N )
                    * ( sumSquaredVoxelValuesReflected - sumVoxelValuesReflected * sumVoxelValuesReflected / N ) ) );
    return cc;
  }

  SImageType::Pointer GetMSPCenteredImage(void)
  {
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

  void Update()
  {
    itk::NumberToString<double> doubleToString;

    if( this->m_DoPowell )
      {
      try
        {
        this->m_Optimizer->StartOptimization();
        }
      catch (itk::ExceptionObject &e)
        {
        std::cout << "Exception thrown ! " << std::endl;
        std::cout << "An error occurred during Optimization" << std::endl;
        std::cout << "Location    = " << e.GetLocation() << std::endl;
        std::cout << "Description = " << e.GetDescription() << std::endl;
        //return EXIT_FAILURE;
        }
      this->m_params = this->m_Optimizer->GetCurrentPosition();
      }
    this->m_cc = this->GetValue();

    std::cout << doubleToString(this->m_params[0] * 180.0 / vnl_math::pi) << " "
              << doubleToString(this->m_params[1] * 180.0 / vnl_math::pi) << " "
              << this->m_params[2] << " cc= "
              << doubleToString(this->m_cc) << " iters= " << this->m_Optimizer->GetCurrentIteration()
              << std::endl;
  }

  itkSetMacro(DoPowell,bool);
  itkGetConstMacro(DoPowell,bool);

  void SetCenterOfHeadMass(const SImageType::PointType & centerOfHeadMass)
  {
    this->m_CenterOfHeadMass = centerOfHeadMass;
    this->m_CenterOfHeadMassIsSet = true;
  }

private:

  const SImageType::PointType & GetCenterOfHeadMass() const
  {
    if( !m_CenterOfHeadMassIsSet )
      {
      itkGenericExceptionMacro(<< "ERROR: m_CenterOfHeadMass is not set!" << std::endl);
      }
    return this->m_CenterOfHeadMass;
  }

  typedef itk::ResampleImageFilter<SImageType, SImageType> ResampleFilterType;

  ParametersType                    m_params;
  SImageType::Pointer               m_OriginalImage;
  SImageType::Pointer               m_ResamplerReferenceImage;
  SImageType::PointType             m_CenterOfHeadMass;
  bool                              m_CenterOfHeadMassIsSet;
  SImageType::PixelType             m_BackgroundValue;
  OptimizerPointer                  m_Optimizer;
  bool                              m_DoPowell;
  LinearInterpolatorType::Pointer   m_imInterp;
  double                            m_cc;
  bool                              m_HasLocalSupport;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkReflectiveCorrelationCenterToImageMetric.hxx"
#endif

#endif
