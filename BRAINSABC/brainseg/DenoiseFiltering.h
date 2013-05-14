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
