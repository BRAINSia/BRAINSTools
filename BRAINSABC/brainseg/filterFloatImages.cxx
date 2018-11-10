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
#include "itkCurvatureFlowImageFilter.h"

#include "itkGradientAnisotropicDiffusionImageFilter.h"

#include "filterFloatImages.h"

void
filterFloatImages(
  std::vector<itk::Image<float, 3>::Pointer> & images,
  std::string & method,
  unsigned int iters,
  double dt)
{
  using FloatImageType = itk::Image<float, 3>;

  if( method.compare("GradientAnisotropicDiffusion") == 0 )
    {
    std::cout << "Gradient Anisotropic Diffusion" << std::endl;
    using AnisoFilterType = itk::GradientAnisotropicDiffusionImageFilter
      <FloatImageType, FloatImageType>;
    for( unsigned i = 0; i < images.size(); i++ )
      {
      AnisoFilterType::Pointer anisofilt = AnisoFilterType::New();

      // if (debugflag)
      //  anisofilt->DebugOn();

      anisofilt->SetNumberOfIterations(iters);
      anisofilt->SetTimeStep(dt);
      anisofilt->SetInput(images[i]);
      anisofilt->Update();
      images[i] = anisofilt->GetOutput();
      }
    }
  // else if? Default is curvature flow
  else
    {
    std::cout << "K flow" << std::endl;
    using CurvatureFilterType = itk::CurvatureFlowImageFilter
      <FloatImageType, FloatImageType>;
    for( unsigned int k = 0; k < images.size(); k++ )
      {
      CurvatureFilterType::Pointer cfilt = CurvatureFilterType::New();

      cfilt->SetNumberOfIterations(iters);
      cfilt->SetTimeStep(dt);
      cfilt->SetInput(images[k]);
      cfilt->Update();

      images[k] = cfilt->GetOutput();
      }
    }
}
