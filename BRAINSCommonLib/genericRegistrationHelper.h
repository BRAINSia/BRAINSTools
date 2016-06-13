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
/*=========================================================================
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile$
 *  Language:  C++
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef __genericRegistrationHelper_h
#define __genericRegistrationHelper_h

#include "BRAINSCommonLib.h"

#include "itkImageToImageMetricv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"

#include "itkImageRegistrationMethodv4.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImage.h"
#include "itkDataObjectDecorator.h"

#include "itkCenteredTransformInitializer.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>

#include "itkMultiThreader.h"
#include "itkResampleImageFilter.h"

#include "itkIntensityWindowingImageFilter.h"

#ifdef USE_DebugImageViewer
#include "DebugImageViewerClient.h"
#include "itkLinearInterpolateImageFunction.h"
#include "Imgmath.h"
#endif

#include "itkIO.h"

#include <cstdio>

enum
  {
  Dimension = 3,
  MaxInputDimension = 4
  };

//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"

namespace BRAINSFit
{
template <class JointPDFType>
void MakeDebugJointHistogram(const std::string & debugOutputDirectory, const typename JointPDFType::Pointer myHistogram,
                             const int globalIteration,
                             const int currentIteration)
{
  std::stringstream fn("");

  fn << debugOutputDirectory << "/DEBUGHistogram_"
     << std::setw( 4 ) << std::setfill( '0' ) << globalIteration
     << "_"
     << std::setw( 4 ) << std::setfill( '0' ) << currentIteration
     << ".png";

  itk::ImageRegionConstIterator<JointPDFType> origIter(myHistogram, myHistogram->GetLargestPossibleRegion() );
  origIter.GoToBegin();
  unsigned int nonZeroCount   = 0;
  float        nonZeroAverage = 0;

  while( !origIter.IsAtEnd() )
    {
    const float currValue = origIter.Get();
    if( currValue > 0 )
      {
      nonZeroAverage += currValue;
      ++nonZeroCount;
      }
    ++origIter;
    }

  nonZeroAverage /= static_cast<float>(nonZeroCount);

  typedef itk::Image<unsigned short, 2> PNGImageType;
  typename PNGImageType::Pointer myOut = PNGImageType::New();
  myOut->CopyInformation(myHistogram);
  myOut->SetRegions(myHistogram->GetLargestPossibleRegion() );
  myOut->Allocate();
  myOut->FillBuffer(0U);
  itk::ImageRegionIterator<PNGImageType> pngIter(myOut, myOut->GetLargestPossibleRegion() );
  pngIter.GoToBegin();
  origIter.GoToBegin();
  while( !pngIter.IsAtEnd() )
    {
    const float MAX_VALUE = std::numeric_limits<unsigned short>::max();
    const float scaleFactor = 0.66 * MAX_VALUE / nonZeroAverage;
    const float currValue = (origIter.Get() ) * scaleFactor;
    if( currValue < 0 )
      {
      pngIter.Set( 0 );
      }
    else if( currValue > MAX_VALUE )
      {
      pngIter.Set( MAX_VALUE );
      }
    else
      {
      pngIter.Set( currValue );
      }
    ++pngIter;
    ++origIter;
    }

  itkUtil::WriteImage<PNGImageType>( myOut, fn.str() );
  std::cout << "Writing jointPDF: " << fn.str() << std::endl;
}

template <typename TOptimizer, typename TTransform, typename TImage>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate Self;
  typedef  itk::Command           Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro(Self);

  typedef          TOptimizer  OptimizerType;
  typedef const OptimizerType *OptimizerPointer;

  typedef const typename OptimizerType::ParametersType OptimizerParametersType;

  typedef typename itk::MattesMutualInformationImageToImageMetricv4<TImage, TImage> MattesMutualInformationMetricType;
  void SetDisplayDeformedImage(bool x)
  {
    m_DisplayDeformedImage = x;
#ifdef USE_DebugImageViewer
    m_DebugImageDisplaySender.SetEnabled(x);
#endif
  }

  void SetPromptUserAfterDisplay(bool x)
  {
    m_PromptUserAfterDisplay = x;
#ifdef USE_DebugImageViewer
    m_DebugImageDisplaySender.SetPromptUser(x);
#endif
  }

  void SetPrintParameters(bool x)
  {
    m_PrintParameters = x;
  }

  void SetFixedImage(typename TImage::Pointer fixed)
  {
    m_FixedImage = fixed;
  }

  void SetMovingImage(typename TImage::Pointer moving)
  {
    m_MovingImage = moving;
  }

  void SetTransform(typename TTransform::Pointer & xfrm)
  {
    m_Transform = xfrm;
  }

  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
  {
    Execute( (const itk::Object *)caller, event );
  }

  typename TImage::Pointer Transform(typename TTransform::Pointer & xfrm)
  {
    typedef typename itk::LinearInterpolateImageFunction<TImage, double> InterpolatorType;
    typename InterpolatorType::Pointer interp = InterpolatorType::New();
    typedef typename itk::ResampleImageFilter<TImage, TImage> ResampleImageFilter;
    typename ResampleImageFilter::Pointer resample = ResampleImageFilter::New();
    resample->SetInput(m_MovingImage);
    resample->SetTransform(xfrm);
    resample->SetInterpolator(interp);
    resample->SetOutputParametersFromImage(m_FixedImage);
    resample->SetOutputDirection( m_FixedImage->GetDirection() );
    resample->SetDefaultPixelValue(0);
    resample->Update();
    typename TImage::Pointer rval = resample->GetOutput();
    return rval;
  }

  void Execute(const itk::Object *object, const itk::EventObject & event) ITK_OVERRIDE
  {
    OptimizerPointer optimizer = dynamic_cast<OptimizerPointer>( object );

    if( optimizer == ITK_NULLPTR )
      {
      itkGenericExceptionMacro("fail to convert to Optimizer Pointer");
      }
    if( !itk::IterationEvent().CheckEvent(&event) )
      {
      return;
      }

    OptimizerParametersType parms =
      optimizer->GetCurrentPosition();
    int  psize = parms.GetNumberOfElements();
    bool parmsNonEmpty = false;
    for( int i = 0; i < psize; ++i )
      {
      if( parms[i] != 0.0 )
        {
        parmsNonEmpty = true;
        break;
        }
      }

    if( m_PrintParameters )
      {
      std::cout << std::setw(  4 ) << std::setfill( ' ' ) << optimizer->GetCurrentIteration() << "   ";
      std::cout << std::setw( 10 ) << std::setfill( ' ' ) << optimizer->GetValue() << "   ";
      if( parmsNonEmpty ) // For BSpline tranform with large parameters space (>15), every (psize/15)th parameter is printed.
        {
        std::cout << "[ ";
        int i = 0;
        while( i<psize )
          {
          std::cout << parms[i] << " ";
          i += std::max(1,(psize/15));
          }
        std::cout << "]" << std::endl;
        }
      }
    //
    // GenerateHistogram
    // TODO: KENT:  BRAINSFit tools need to define a common output directory for
    // all debug images to be written.
    //             by default the it should be the same as the outputImage, and
    // if that does not exists, then it
    //             should default to the same directory as the outputTransform,
    // or it should be specified by the
    //             user on the command line.
    //             The following function should only be called when BRAINSFit
    // command line tool is called with
    //             --debugLevel 7 or greater, and it should write the 3D
    // JointPDF image to the debugOutputDirectory
    //             location.
    std::string debugOutputDirectory("");
    if( optimizer->GetCurrentIteration() < 5 // Only do first 4 of each iteration
        &&
        itksys::SystemTools::GetEnv("DEBUG_JOINTHISTOGRAM_DIR", debugOutputDirectory) && debugOutputDirectory != "" )
      {
      // Special BUG work around for MMI metric
      // that does not work in multi-threaded mode
      // typedef itk::MattesMutualInformationImageToImageMetricv4<TImage,TImage> MattesMutualInformationMetricType;
      static int TransformIterationCounter = 0;
      typename MattesMutualInformationMetricType::Pointer test_MMICostMetric =
        dynamic_cast<MattesMutualInformationMetricType *>(this->m_ObserverCostMetricObject.GetPointer() );
      if( test_MMICostMetric.IsNotNull() )
        {
        typedef typename MattesMutualInformationMetricType::PDFValueType PDFValueType;
        typedef itk::Image<PDFValueType, 2>                              JointPDFType;
        const typename JointPDFType::Pointer myHistogram = test_MMICostMetric->GetJointPDF();
        MakeDebugJointHistogram<JointPDFType>(debugOutputDirectory, myHistogram, TransformIterationCounter,
                                              optimizer->GetCurrentIteration() );
        if( optimizer->GetCurrentIteration()  == 0 )
          {
          TransformIterationCounter += 10000;
          }
        }
      }

#ifdef USE_DebugImageViewer
    if( m_DisplayDeformedImage )
      {
      if( parmsNonEmpty )
        {
        m_Transform->SetParametersByValue(parms);
        }
      // else, if it is a vnl optimizer wrapper, i.e., the BSpline optimizer,
      // the only hint you get
      // is in the transform object used by the optimizer, so don't erase it,
      // use it.
      typename TImage::Pointer transformResult =
        this->Transform(m_Transform);
      m_DebugImageDisplaySender.SendImage<TImage>(transformResult, 1);
      typename TImage::Pointer diff = Isub<TImage>(m_FixedImage, transformResult);

      m_DebugImageDisplaySender.SendImage<TImage>(diff, 2);
      }
#endif
  }

  void SetMetricObject(typename MattesMutualInformationMetricType::Pointer metric_Object)
  {
    // NOTE:  Returns NULL if not MattesMutualInformationImageToImageMetric
    this->m_ObserverCostMetricObject = dynamic_cast<MattesMutualInformationMetricType *>(metric_Object.GetPointer() );
    // NO NEED FOR CHECKING IF DYNAMIC CAST WORKED, invalid cast is OK with a
    // NULL return
  }

protected:
  CommandIterationUpdate() : m_DisplayDeformedImage(false),
    m_PromptUserAfterDisplay(false),
    m_PrintParameters(true),
    m_MovingImage(ITK_NULLPTR),
    m_FixedImage(ITK_NULLPTR),
    m_Transform(ITK_NULLPTR),
    m_ObserverCostMetricObject(ITK_NULLPTR)
  {
  }

private:
  bool m_DisplayDeformedImage;
  bool m_PromptUserAfterDisplay;
  bool m_PrintParameters;

  typename TImage::Pointer m_MovingImage;
  typename TImage::Pointer m_FixedImage;
  typename TTransform::Pointer m_Transform;

  typename MattesMutualInformationMetricType::Pointer m_ObserverCostMetricObject;

#ifdef USE_DebugImageViewer
  DebugImageViewerClient m_DebugImageDisplaySender;
#endif
};
} // end namespace BRAINSFit

namespace itk
{
// TODO:  Remove default MetricType here, and force a choice
template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
class MultiModal3DMutualRegistrationHelper : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef MultiModal3DMutualRegistrationHelper Self;
  typedef ProcessObject                        Superclass;
  typedef SmartPointer<Self>                   Pointer;
  typedef SmartPointer<const Self>             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiModal3DMutualRegistrationHelper, ProcessObject);

  typedef          TFixedImage                  FixedImageType;
  typedef typename FixedImageType::ConstPointer FixedImageConstPointer;
  typedef typename FixedImageType::Pointer      FixedImagePointer;

  typedef          TMovingImage                  MovingImageType;
  typedef typename MovingImageType::ConstPointer MovingImageConstPointer;
  typedef typename MovingImageType::Pointer      MovingImagePointer;

  typedef          TTransformType         TransformType;
  typedef typename TransformType::Pointer TransformPointer;

  /** Type for the output: Using Decorator pattern for enabling
    *  the Transform to be passed in the data pipeline */
  typedef DataObjectDecorator<TransformType>    TransformOutputType;
  typedef typename TransformOutputType::Pointer TransformOutputPointer;
  typedef typename TransformOutputType::ConstPointer
    TransformOutputConstPointer;

  /** Constants for the image dimensions */
  itkStaticConstMacro(FixedImageDimension, unsigned int, FixedImageType::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, unsigned int, MovingImageType::ImageDimension);

  typedef typename itk::ObjectToObjectMultiMetricv4< FixedImageDimension,
                                                     MovingImageDimension,
                                                     FixedImageType,
                                                     double>                  MultiMetricType;
  typedef typename itk::ImageToImageMetricv4< FixedImageType,
                                              MovingImageType,
                                              FixedImageType,
                                              double >                        ImageMetricType;

  typedef itk::CompositeTransform<double, MovingImageDimension>     CompositeTransformType;

  typedef          TOptimizer                    OptimizerType;
  typedef const OptimizerType *                  OptimizerPointer;
  typedef typename OptimizerType::ScalesType     OptimizerScalesType;
  typedef typename OptimizerType::ParametersType OptimizerParametersType;

  typedef itk::GradientDescentOptimizerBasev4Template< double > GenericOptimizerType;

  typedef ImageRegistrationMethodv4< FixedImageType,
                                     MovingImageType >       RegistrationType;
  typedef typename RegistrationType::Pointer                 RegistrationPointer;

  typedef itk::AffineTransform<double, 3>                                  AffineTransformType;
  typedef itk::ImageRegistrationMethodv4< FixedImageType, MovingImageType> AffineRegistrationType;
  typedef typename AffineRegistrationType::MetricSamplingStrategyType      SamplingStrategyType;

  typedef itk::CenteredTransformInitializer<
      TransformType,
      FixedImageType,
      MovingImageType>    TransformInitializerType;

  typedef itk::ResampleImageFilter<
      MovingImageType,
      FixedImageType>     ResampleFilterType;

  /** Initialize by setting the interconnects between the components. */
  virtual void Initialize(void); // throw ( ExceptionObject );

  /** Method that initiates the registration. */
  void Update(void) ITK_OVERRIDE;

  /** Set/Get the Fixed image. */
  void SetFixedImage(FixedImagePointer fixedImage);

  itkGetConstObjectMacro(FixedImage, FixedImageType);

  /** Set/Get the Fixed image2. */
  void SetFixedImage2(FixedImagePointer fixedImage2);

  itkGetConstObjectMacro(FixedImage2, FixedImageType);

  /** Set/Get the Moving image. */
  void SetMovingImage(MovingImagePointer movingImage);

  itkGetConstObjectMacro(MovingImage, MovingImageType);

  /** Set/Get the Moving image2. */
  void SetMovingImage2(MovingImagePointer movingImage2);

  itkGetConstObjectMacro(MovingImage2, MovingImageType);

  /** Set/Get the InitialTransfrom. */
  void SetInitialTransform(typename CompositeTransformType::Pointer initialTransform);

  /** Set/Get the Transfrom. */
  itkSetObjectMacro(Transform, TransformType);
  typename CompositeTransformType::Pointer GetTransform(void);

  itkSetObjectMacro(CostMetricObject, MetricType);
  itkGetConstObjectMacro(CostMetricObject, MetricType);

  itkSetMacro(SamplingPercentage,            double);
  itkSetMacro(NumberOfHistogramBins,         unsigned int);
  itkSetMacro(NumberOfIterations,            unsigned int);
  itkSetMacro(RelaxationFactor,              double);
  itkSetMacro(MaximumStepLength,             double);
  itkSetMacro(MinimumStepLength,             double);
  itkSetMacro(TranslationScale,              double);
  itkSetMacro(ReproportionScale,             double);
  itkSetMacro(SkewScale,                     double);
  itkSetMacro(BackgroundFillValue,           double);
  itkSetMacro(DisplayDeformedImage,          bool);
  itkSetMacro(PromptUserAfterDisplay,        bool);
  itkGetConstMacro(FinalMetricValue,         double);
  itkGetConstMacro(ActualNumberOfIterations, unsigned int);
  itkSetMacro(ObserveIterations,        bool);
  itkGetConstMacro(ObserveIterations,        bool);

  itkSetMacro(SamplingStrategy,SamplingStrategyType);
  itkGetConstMacro(SamplingStrategy,SamplingStrategyType);

  /** Returns the transform resulting from the registration process  */
  const TransformOutputType * GetOutput() const;

  /** Make a DataObject of the correct type to be used as the specified
    * output. */
  using Superclass::MakeOutput;
  virtual DataObjectPointer MakeOutput(unsigned int idx);

  /** Method to return the latest modified time of this object or
    * any of its cached ivars */
  unsigned long GetMTime() const ITK_OVERRIDE;

protected:
  MultiModal3DMutualRegistrationHelper();
  virtual ~MultiModal3DMutualRegistrationHelper()
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  /** Method invoked by the pipeline in order to trigger the computation of
    * the registration. */
  void  GenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiModal3DMutualRegistrationHelper);

  FixedImagePointer  m_FixedImage;
  MovingImagePointer m_MovingImage;
  FixedImagePointer  m_FixedImage2;
  MovingImagePointer m_MovingImage2;

  typename CompositeTransformType::Pointer m_CompositeTransform;
  TransformPointer                         m_Transform;
  //
  // make sure parameters persist until after they're used by the transform

  RegistrationPointer m_Registration;

  typename MetricType::Pointer  m_CostMetricObject;

  double       m_SamplingPercentage;
  unsigned int m_NumberOfHistogramBins;
  unsigned int m_NumberOfIterations;
  double       m_RelaxationFactor;
  double       m_MaximumStepLength;
  double       m_MinimumStepLength;
  double       m_TranslationScale;
  double       m_ReproportionScale;
  double       m_SkewScale;
  double       m_BackgroundFillValue;
  unsigned int m_ActualNumberOfIterations;
  bool         m_DisplayDeformedImage;
  bool         m_PromptUserAfterDisplay;
  double       m_FinalMetricValue;
  bool         m_ObserveIterations;

  SamplingStrategyType m_SamplingStrategy;

  ModifiedTimeType m_InternalTransformTime;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "genericRegistrationHelper.hxx"
#endif

#endif
