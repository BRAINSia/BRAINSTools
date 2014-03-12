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
#ifndef __DenoiseFiltering_h
#define __DenoiseFiltering_h

#include "itkCurvatureFlowImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"

// TODO:  BRAINSFit currently has an option for doing a median filter to remove
// noise,
//       this should be expanded to allow for a richer set of pre-filter
// denoising operations
//       such as the anisotropic diffusion denoising options.
//       This should be made into a strategy pattern class that is dynamically
// selectable at runtime.
//

template <class TInputImageType>
// TODO:  Input and outputs should be templated separately?
typename TInputImageType::Pointer DenoiseFiltering(
  typename TInputImageType::Pointer img,
  const std::string & PrefilteringMethod,    // Select the type of denoising to do
  const unsigned int PrefilteringIterations, // Only used in AD and CF
  const double PrefilteringTimeStep,         // Only used in AD and CF
  const std::vector<unsigned int> &          // gridSize   //Only used in median filtering
  )
{
  typename TInputImageType::Pointer denoisedImage = NULL;
  if( PrefilteringMethod.compare("GradientAnisotropicDiffusion") == 0 && PrefilteringIterations >  0 )
    {
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    std::cout << "Prefiltering with " << PrefilteringMethod << " (Iters=" << PrefilteringIterations
              << ",TimeStep=" << PrefilteringTimeStep << ")  gridSize=(" << ")" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    typedef itk::GradientAnisotropicDiffusionImageFilter<TInputImageType, TInputImageType> AnisoFilterType;
    typename AnisoFilterType::Pointer anisofilt = AnisoFilterType::New();

    anisofilt->SetInput(img);
    anisofilt->SetConductanceParameter(1);
    anisofilt->SetNumberOfIterations(PrefilteringIterations);
    anisofilt->SetTimeStep(PrefilteringTimeStep);
    anisofilt->Update();

    denoisedImage = anisofilt->GetOutput();
    }
  else if( PrefilteringMethod.compare("CurvatureFlow") == 0   && PrefilteringIterations > 0 )
    {
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    std::cout << "Prefiltering with " << PrefilteringMethod << " (Iters=" << PrefilteringIterations
              << ",TimeStep=" << PrefilteringTimeStep << ")  gridSize=(" << ")" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    typedef itk::CurvatureFlowImageFilter<TInputImageType, TInputImageType> CurvatureFilterType;
    typename CurvatureFilterType::Pointer cfilt = CurvatureFilterType::New();
    cfilt->SetInput(img);
    cfilt->SetNumberOfIterations(PrefilteringIterations);
    cfilt->SetTimeStep(PrefilteringTimeStep);
    cfilt->Update();

    denoisedImage = cfilt->GetOutput();
    }
  else if( PrefilteringMethod.compare("MedianFilter") == 0 )
    {
    // TODO:  Kent put a median filter in here, with filter radius equal to 1.
    }

  if( denoisedImage.IsNotNull() )
    {
    return denoisedImage;
    }
  // else No filtering.
  return img;
}


#endif
