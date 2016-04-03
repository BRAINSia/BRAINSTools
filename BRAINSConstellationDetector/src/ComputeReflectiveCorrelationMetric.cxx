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
// Author: Ali Ghayoor

#include <iostream>
#include "itkIO.h"
#include "itkImageFileReader.h"
#include "itkMultiResolutionPyramidImageFilter.h"

#include "landmarksConstellationCommon.h"
#include "StandardizeMaskIntensity.h"
#include "itkFindCenterOfBrainFilter.h"

#define WRITE_CSV_FILE
#include "itkReflectiveCorrelationCenterToImageMetric.h"
#undef WRITE_CSV_FILE

#include "ComputeReflectiveCorrelationMetricCLP.h"


typedef itk::FindCenterOfBrainFilter<SImageType>                        FindCenterFilter;
typedef Rigid3DCenterReflectorFunctor< itk::PowellOptimizerv4<double> > ReflectionFunctorType;
typedef ReflectionFunctorType::ParametersType                           ParametersType;
typedef itk::CastImageFilter<DImageType3D, SImageType>                  CasterType;

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  // load image
  std::cout << "\nLoading image..." << std::endl;
  // Input image is read as a double image;
  // then it is rescaled to a specific dynamic range;
  // Finally it is cast to a Short type image.
  typedef itk::ImageFileReader<DImageType3D> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << " Error while reading image file( s ) with ITK:\n "
              << err << std::endl;
    }

  DImageType3D::Pointer rescaledInputVolume =
    StandardizeMaskIntensity<DImageType3D, ByteImageType>(reader->GetOutput(),
                                                          ITK_NULLPTR,
                                                          0.0005, 1.0 - 0.0005,
                                                          1, 0.95 * MAX_IMAGE_OUTPUT_VALUE,
                                                          0, MAX_IMAGE_OUTPUT_VALUE);

  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( rescaledInputVolume );
  caster->Update();
  SImageType::Pointer originalImage = caster->GetOutput();

#if 0 //HACK removed for speed testing.
  // Find center of head mass
  std::cout << "\nFinding center of head mass..." << std::endl;
  FindCenterFilter::Pointer findCenterFilter = FindCenterFilter::New();
  findCenterFilter->SetInput( originalImage );
  findCenterFilter->SetAxis( 2 );
  findCenterFilter->SetOtsuPercentileThreshold( 0.01 );
  findCenterFilter->SetClosingSize( 7 );
  findCenterFilter->SetHeadSizeLimit( 700 );
  findCenterFilter->SetBackgroundValue( 0 );
  findCenterFilter->Update();
  SImagePointType centerOfHeadMass = findCenterFilter->GetCenterOfBrain();
#else
  SImagePointType centerOfHeadMass;
  centerOfHeadMass[0] = -1.7824640167960075;
  centerOfHeadMass[1] = -11.914627390460533;
  centerOfHeadMass[2] = 32.393109592514605;
#endif

  PyramidFilterType::Pointer MyPyramid = MakeOneLevelPyramid( originalImage );
  SImageType::Pointer inputImage = MyPyramid->GetOutput(0); // one-eighth image

  ReflectionFunctorType::Pointer reflectionFunctor = ReflectionFunctorType::New();
  reflectionFunctor->InitializeImage( originalImage ); // Initialize from the high-res image
  reflectionFunctor->SetDownSampledReferenceImage(inputImage);
  reflectionFunctor->SetCenterOfHeadMass(centerOfHeadMass); // center of head mass must be passed

  // optimal parameters
  ParametersType opt_params;
  opt_params.set_size(ReflectionFunctorType::SpaceDimension);
  opt_params.fill(0.0);
  reflectionFunctor->SetParameters(opt_params);
  reflectionFunctor->SetDoPowell(false);
  reflectionFunctor->Update();
  double opt_cc = reflectionFunctor->GetValue();

#if 0
  const double HA_range = 45.0;
  const double BA_range = 45.0;
  const double LR_range = 5.0;

  const double HA_stepsize = 5.0; // degree
  const double BA_stepsize = 5.0; // degree
  const double LR_stepsize = 1.0; // mm
#else
  const double HA_range = 10.0;
  const double BA_range = 25.0;
  const double LR_range = 5.0;

  const double HA_stepsize = 5.0; // degree
  const double BA_stepsize = 5.0; // degree
  const double LR_stepsize = 2.0; // mm
#endif


  reflectionFunctor->DoExhaustiveSearch(opt_params, opt_cc,
                                        HA_range, BA_range, LR_range,
                                        HA_stepsize, BA_stepsize, LR_stepsize,
                                        outputCSVFile);

  std::cout << "Optimize parameters by exhaustive search: [" << opt_params[0] << "," << opt_params[1] << "," << opt_params[2] << "]" << std::endl;
  std::cout << "Optimize metric value by exhaustive search: " << opt_cc << std::endl;

/*
  // Now compare find the optimal parameters using Powell Optimizer
  ReflectionFunctorType::Pointer reflectionFunctor2 = ReflectionFunctorType::New();
  reflectionFunctor2->SetDownSampledReferenceImage(inputImage);
  reflectionFunctor2->Initialize();
  reflectionFunctor2->Update();
  ParametersType powell_params = reflectionFunctor2->GetParameters();
  double powell_cc = reflectionFunctor2->GetValue();

  std::cout << "Optimize parameters by Powell search: [" << powell_params[0] << "," << powell_params[1] << "," << powell_params[2] << "]" << std::endl;
  std::cout << "Optimize metric value by Powell search: " << powell_cc << std::endl;
*/
  // here compare opt_params with input baseline params to return failure or success.

  return EXIT_SUCCESS;
}
