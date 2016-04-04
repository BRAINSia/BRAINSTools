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


#include "itkTimeProbe.h"

#define WRITE_CSV_FILE
#include "itkReflectiveCorrelationCenterToImageMetric.h"
#undef WRITE_CSV_FILE

#include "ComputeReflectiveCorrelationMetricCLP.h"



int main( int argc, char * argv[] ) {
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  // load image
  std::cout << "\nLoading image..." << std::endl;
  // Input image is read as a double image;
  // then it is rescaled to a specific dynamic range;
  // Finally it is cast to a Short type image.
  typedef itk::ImageFileReader<DImageType3D> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputVolume);
  try {
    reader->Update();
  }
  catch (itk::ExceptionObject &err) {
    std::cerr << " Error while reading image file( s ) with ITK:\n "
    << err << std::endl;
  }

  DImageType3D::Pointer rescaledInputVolume =
      StandardizeMaskIntensity<DImageType3D, ByteImageType>(reader->GetOutput(),
                                                            ITK_NULLPTR,
                                                            0.0005, 1.0 - 0.0005,
                                                            1, 0.95 * MAX_IMAGE_OUTPUT_VALUE,
                                                            0, MAX_IMAGE_OUTPUT_VALUE);

  typedef itk::CastImageFilter<DImageType3D, SImageType>                  CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(rescaledInputVolume);
  caster->Update();
  SImageType::Pointer originalImage = caster->GetOutput();

  PyramidFilterType::Pointer MyPyramid = MakeOneLevelPyramid(originalImage);
  SImageType::Pointer inputImage = MyPyramid->GetOutput(0); // one-eighth image

  typedef Rigid3DCenterReflectorFunctor< itk::PowellOptimizerv4<double> > ReflectionFunctorType;
  typedef ReflectionFunctorType::ParametersType                           ParametersType;

  ReflectionFunctorType::Pointer reflectionFunctor = ReflectionFunctorType::New();
  reflectionFunctor->InitializeImage(inputImage);

  // optimal parameters
  ParametersType opt_params;
  opt_params.set_size(ReflectionFunctorType::SpaceDimension);
  opt_params.fill(0.0);
  reflectionFunctor->SetParameters(opt_params);
  reflectionFunctor->SetDoPowell(false);
  reflectionFunctor->Update();
  double opt_cc = reflectionFunctor->GetValue();


  std::vector<std::string> suffix(4);
  suffix[0]="0";
  suffix[1]="1";
  suffix[2]="2";

  std::vector<double> Angle_Range(3);
  Angle_Range[0]=45.0;
  Angle_Range[1]=2.5;
  Angle_Range[2]=0.5;

  std::vector<double> Angle_Stepsizes(3);
  Angle_Stepsizes[0] = 5.0;
  Angle_Stepsizes[1] = 0.5;
  Angle_Stepsizes[2] = 0.25;

  std::vector<double> Offset_Range(3);
  Offset_Range[0]=15.0;
  Offset_Range[1]=1.5;
  Offset_Range[2]=0.25;

  std::vector<double> Offset_Stepsizes(3);
  Offset_Stepsizes[0] = 3.0;
  Offset_Stepsizes[1] = 0.5;
  Offset_Stepsizes[2] = .25;

  for (unsigned int resolutionIter = 0; resolutionIter <= 2; ++resolutionIter )
  {
    const double HA_range =  Angle_Range[resolutionIter];
    const double BA_range =  Angle_Range[resolutionIter];
    const double LR_range =  Offset_Range[resolutionIter];

    const double HA_stepsize = Angle_Stepsizes[resolutionIter]; // degree
    const double BA_stepsize = Angle_Stepsizes[resolutionIter]; // degree
    const double LR_stepsize = Offset_Stepsizes[resolutionIter]; // mm

    std::cout << "RANGE: " << HA_range << " at " << HA_stepsize << std::endl;
    std::cout << "LR: " << LR_range << " at " << LR_stepsize << std::endl;
    itk::TimeProbe clock;
    clock.Start();
    reflectionFunctor->DoExhaustiveSearch(opt_params, opt_cc,
                                        HA_range, BA_range, LR_range,
                                        HA_stepsize, BA_stepsize, LR_stepsize,
                                        outputCSVFile+suffix[resolutionIter]);
    clock.Stop();
    std::cout << "Mean: " << clock.GetMean() << std::endl;
    std::cout << "Total: " << clock.GetTotal() << std::endl;

    const double degree_to_rad = vnl_math::pi / 180.0;

    std::cout << "Optimize parameters by exhaustive search: [" << opt_params[0]/degree_to_rad << "," <<
      opt_params[1]/degree_to_rad << "," << opt_params[2] << "]" << std::endl;
    std::cout << "Optimize metric value by exhaustive search: " << opt_cc << std::endl;
  }

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
