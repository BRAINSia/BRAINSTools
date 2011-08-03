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
  const std::string PrefilteringMethod,      // Select the type of denoising to
                                             // do
  const unsigned int PrefilteringIterations, // Only used in AD and CF
  const double PrefilteringTimeStep,         // Only used in AD and CF
  const std::vector<unsigned int>            // gridSize   //Only used in median
                                             // filtering
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

#if 0
< !--Standard options for the denoising phase of a program-- >
< integer - vector >
< name > medianFilterSize</ name>
< longflag > medianFilterSize</ longflag>
< label > Median Filter Size</ label>
< description > The radius for the optional MedianImageFilter preprocessing in all 3 directions.< / description >
< default > 0, 0, 0 < / default >
< / integer - vector >
< integer >
< name > filterIteration</ name>
< description > Filter iterations</ description>
< label > Filter Iterations</ label>
< longflag > filterIteration</ longflag>
< default > 10 < / default >
< constraints >
< minimum > 0 < / minimum >
< maximum > 50 < / maximum >
< step > 1 < / step >
< / constraints >
< / integer >
< float >
< name > filterTimeStep</ name>
< description > Filter time step</ description>
< label > Filter Time Step</ label>
< longflag > filterTimeStep</ longflag>
< default > 0.01 < / default >
< constraints >
< minimum > 0 < / minimum >
< maximum > 0.5 < / maximum >
< step > 0.01 < / step >
< / constraints >
< / float >
< string - enumeration >
< name > filterMethod</ name>
< label > Filter Method</ label>
< longflag > filterMethod</ longflag>
< description > Filter method for preprocessing of registration</ description>
< element > None</ element>
< element > CurvatureFlow</ element>
< element > GradientAnisotropicDiffusion</ element>
< element > Median</ element>
< default > None</ default>
< / string - enumeration >
#endif

#endif
