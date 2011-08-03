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
  typedef itk::Image<float, 3> FloatImageType;

  if( method.compare("GradientAnisotropicDiffusion") == 0 )
    {
    std::cout << "Gradient Anisotropic Diffusion" << std::endl;
    typedef itk::GradientAnisotropicDiffusionImageFilter
      <FloatImageType, FloatImageType>
      AnisoFilterType;
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
    typedef itk::CurvatureFlowImageFilter
      <FloatImageType, FloatImageType> CurvatureFilterType;
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
