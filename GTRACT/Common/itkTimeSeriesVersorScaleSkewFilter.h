/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkTimeSeriesVersorScaleSkewFilter_h
#define __itkTimeSeriesVersorScaleSkewFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include <itkExtractImageFilter.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkResampleImageFilter.h>
#include <itkImageRegistrationMethod.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkScaleSkewVersor3DTransform.h>
#include <itkCenteredTransformInitializer.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkTransformFactory.h>
#include "gtractCommonWin32.h"

#include <map>
#include <string>

namespace itk
{
/** \class Orient4dImageFilter
 * \brief Permute axes and then flip images as needed to obtain
 *  agreement in coordinateOrientation codes.
 *
 * This class satisfies performs the following steps:
 *    For i in 4th Dimension
 *      ExtractVolume with Extract Image Filter
 *      Orient 3D extracted volume
 *    End
 *
 * It is build upon the ExtractImageFilter and the OrientImageFilter
 */

class GTRACT_COMMON_EXPORT TimeSeriesVersorScaleSkewFilter : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef TimeSeriesVersorScaleSkewFilter Self;
  typedef itk::Object                     Superclass;
  typedef SmartPointer<Self>              Pointer;
  typedef SmartPointer<const Self>        ConstPointer;

  /** Some convenient typedefs. */
  typedef itk::Image<signed short, 4>   InputImageType;
  typedef itk::Image<signed short, 4>   OutputImageType;
  typedef InputImageType::Pointer       InputImagePointer;
  typedef InputImageType::ConstPointer  InputImageConstPointer;
  typedef InputImageType::RegionType    InputImageRegionType;
  typedef InputImageType::SizeType      InputImageSizeType;
  typedef InputImageType::SpacingType   InputImageSpacingType;
  typedef InputImageType::PointType     InputImagePointType;
  typedef InputImageType::PixelType     InputImagePixelType;
  typedef InputImageType::DirectionType InputImageDirectionType;
  typedef InputImageType::IndexType     InputImageIndexType;

  typedef OutputImageType::Pointer      OutputImagePointer;
  typedef OutputImageType::ConstPointer OutputImageConstPointer;
  typedef OutputImageType::RegionType   OutputImageRegionType;
  typedef OutputImageType::PixelType    OutputImagePixelType;
  typedef OutputImageType::IndexType    OutputImageIndexType;

  typedef itk::Image<InputImagePixelType, 3> ExtractImageType;
  typedef ExtractImageType::Pointer          ExtractImagePointer;
  typedef itk::ExtractImageFilter<InputImageType, ExtractImageType>
    ExtractFilterType;
  typedef ExtractFilterType::Pointer  ExtractFilterTypePointer;
  typedef ExtractImageType::IndexType ExtractImageIndexType;

  typedef itk::ScaleSkewVersor3DTransform<double>  TransformType;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::MattesMutualInformationImageToImageMetric<
      ExtractImageType,
      ExtractImageType>        MetricType;

  typedef itk::LinearInterpolateImageFunction<
      ExtractImageType,
      double>         InterpolatorType;

  typedef itk::ImageRegistrationMethod<
      ExtractImageType,
      ExtractImageType>        RegistrationType;

  typedef itk::CenteredTransformInitializer<TransformType,
                                            ExtractImageType,
                                            ExtractImageType
                                            >  TransformInitializerType;
  typedef itk::ResampleImageFilter<
      ExtractImageType,
      ExtractImageType>    ResampleFilterType;

  typedef TransformType::Pointer            TransformTypePointer;
  typedef TransformType::VersorType         VersorType;
  typedef VersorType::VectorType            VectorType;
  typedef MetricType::Pointer               MetricTypePointer;
  typedef OptimizerType::Pointer            OptimizerTypePointer;
  typedef OptimizerType::ParametersType     OptimizerParameterType;
  typedef OptimizerType::ScalesType         OptimizerScalesType;
  typedef InterpolatorType::Pointer         InterpolatorTypePointer;
  typedef RegistrationType::Pointer         RegistrationTypePointer;
  typedef TransformInitializerType::Pointer TransformInitializerTypePointer;
  typedef ResampleFilterType::Pointer       ResampleFilterTypePointer;

#if 0
  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** The dimensions of the input image must equal those of the
      output image. */
  itkConceptMacro( SameDimension,
                   ( Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),
                                            itkGetStaticConstMacro(OutputImageDimension)> ) );

/** The dimension of the input image must be 4. */
  itkConceptMacro( DimensionShouldBe4,
                   ( Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension), 4> ) );
#endif
/** Standard New method. */
  itkNewMacro(Self);

/** Runtime information support. */
  itkTypeMacro(Orient4dImageFilter, itk::Object);

/* SetInput and GetOutput Macros */
  itkSetObjectMacro(Input, InputImageType);
  itkGetConstObjectMacro(Output, OutputImageType);

  itkSetMacro(NumberOfSpatialSamples, int);
  itkSetMacro(NumberOfIterations, int);
  itkSetMacro(TranslationScale, float);
  itkSetMacro(ScalingScale, float);
  itkSetMacro(SkewScale, float);
  itkSetMacro(MaximumStepLength, float);
  itkSetMacro(MinimumStepLength, float);
  itkSetMacro(RelaxationFactor, float);
  itkSetMacro(NumberOfHistogramBins, int);
  itkSetMacro(BaseImage, int);
  itkSetMacro(ResultFile, std::string);

  void Update();

protected:
  TimeSeriesVersorScaleSkewFilter();
  ~TimeSeriesVersorScaleSkewFilter()
  {
  }

private:
  TimeSeriesVersorScaleSkewFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                  // purposely not implemented

// Input and Output Image
  InputImagePointer  m_Input;
  OutputImagePointer m_Output;

// Optional Flip - Used to fix problems with the Direction Cosines
  float                     m_TranslationScale;
  float                     m_ScalingScale;
  float                     m_SkewScale;
  float                     m_MaximumStepLength;
  float                     m_MinimumStepLength;
  float                     m_RelaxationFactor;
  int                       m_NumberOfSpatialSamples;
  int                       m_NumberOfIterations;
  int                       m_NumberOfHistogramBins;
  int                       m_BaseImage;
  std::string               m_ResultFile;
  static const unsigned int m_transformDimension = 3;
};  // end of class
} // end namespace itk

#endif
