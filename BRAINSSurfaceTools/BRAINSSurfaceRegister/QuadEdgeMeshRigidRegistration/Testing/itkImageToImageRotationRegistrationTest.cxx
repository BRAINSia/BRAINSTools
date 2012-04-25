/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageRegistration8.cxx,v $
  Language:  C++
  Date:      $Date: 2008-03-10 19:46:31 $
  Version:   $Revision: 1.37 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImage.h"
#include "itkCenteredTransformInitializer.h"

#define VERSORTRANSFORM
// #undef VERSORTRANSFORM

#ifdef VERSORTRANSFORM
#include "itkVersorTransform.h"
#include "itkVersorTransformOptimizer.h"
#else
#include "itkVersorRigid3DTransform.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#endif

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate Self;
  typedef  itk::Command           Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate()
  {
  };
public:
#ifdef VERSORTRANSFORM
  typedef itk::VersorTransformOptimizer OptimizerType;
#else
  typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
#endif
  typedef   const OptimizerType * OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer =
      dynamic_cast<OptimizerPointer>( object );

    if( !itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }
};

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile ";
    std::cerr << " outputImagefile intialRotationAngle" << std::endl;
    return EXIT_FAILURE;
    }

  const    unsigned int Dimension = 3;
  typedef  float PixelType;

  typedef itk::Image<PixelType, Dimension> FixedImageType;
  typedef itk::Image<PixelType, Dimension> MovingImageType;

#ifdef VERSORTRANSFORM
  typedef itk::VersorTransform<double>  TransformType;
  typedef itk::VersorTransformOptimizer OptimizerType;
#else
  typedef itk::VersorRigid3DTransform<double>  TransformType;
  typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
#endif

  typedef itk::MeanSquaresImageToImageMetric<
      FixedImageType,
      MovingImageType>    MetricType;

  typedef itk::LinearInterpolateImageFunction<
      MovingImageType,
      double>    InterpolatorType;

  typedef itk::ImageRegistrationMethod<
      FixedImageType,
      MovingImageType>    RegistrationType;

  MetricType::Pointer       metric        = MetricType::New();
  OptimizerType::Pointer    optimizer     = OptimizerType::New();
  InterpolatorType::Pointer interpolator  = InterpolatorType::New();
  RegistrationType::Pointer registration  = RegistrationType::New();

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );

  TransformType::Pointer transform = TransformType::New();
  registration->SetTransform( transform );

  typedef itk::ImageFileReader<FixedImageType>  FixedImageReaderType;
  typedef itk::ImageFileReader<MovingImageType> MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );

  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );
  fixedImageReader->Update();

  registration->SetFixedImageRegion(
    fixedImageReader->GetOutput()->GetBufferedRegion() );

  typedef itk::CenteredTransformInitializer<TransformType,
                                            FixedImageType,
                                            MovingImageType
                                            >  TransformInitializerType;

  TransformInitializerType::Pointer initializer =
    TransformInitializerType::New();

  initializer->SetTransform(   transform );
  initializer->SetFixedImage(  fixedImageReader->GetOutput() );
  initializer->SetMovingImage( movingImageReader->GetOutput() );

  initializer->GeometryOn();

  initializer->InitializeTransform();

  typedef TransformType::VersorType VersorType;
  typedef VersorType::VectorType    VectorType;

  VersorType rotation;
  VectorType axis;

  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;

  const double angle = atof( argv[4] );

  rotation.Set(  axis, angle  );

  transform->SetRotation( rotation );

  registration->SetInitialTransformParameters( transform->GetParameters() );

  typedef OptimizerType::ScalesType OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;

#ifndef VERSORTRANSFORM
  optimizerScales[3] = 1.0;
  optimizerScales[4] = 1.0;
  optimizerScales[5] = 1.0;
#endif

  optimizer->SetScales( optimizerScales );

  optimizer->SetMaximumStepLength( 0.2000  );
  optimizer->SetMinimumStepLength( 0.0001 );

  optimizer->SetNumberOfIterations( 200 );

  // Create the Command observer and register it with the optimizer.
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  try
    {
    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  OptimizerType::ParametersType finalParameters =
    registration->GetLastTransformParameters();

  const double versorX  = finalParameters[0];
  const double versorY  = finalParameters[1];
  const double versorZ  = finalParameters[2];

  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();

  const double bestValue = optimizer->GetValue();

  // Print out results
  std::cout << std::endl << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " versor X      = " << versorX  << std::endl;
  std::cout << " versor Y      = " << versorY  << std::endl;
  std::cout << " versor Z      = " << versorZ  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;

  transform->SetParameters( finalParameters );

  TransformType::MatrixType matrix = transform->GetRotationMatrix();
  TransformType::OffsetType offset = transform->GetOffset();

  std::cout << "Matrix = " << std::endl << matrix << std::endl;
  std::cout << "Offset = " << std::endl << offset << std::endl;

  return EXIT_SUCCESS;
}
