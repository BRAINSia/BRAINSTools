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

#ifndef __itkVectorImageRegisterVersorScaleSkewFilter_h
#define __itkVectorImageRegisterVersorScaleSkewFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkVectorImage.h"
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
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkTransformFileWriter.h>
#include <itkCastImageFilter.h>
#include "gtractCommonWin32.h"

#include <map>
#include <string>

namespace itk
{
/** \class VectorImageRegisterVersorScaleSkewFilter
 * \brief Register each index of a vector image to an image
 *        specified by the user. The fixed image is not a vector
 *        image. The result is a new vector image.
 *
 * This class was designed to handle image co-registeration for
 * time series that may be loaded as vector images. This
 * may the case for DWI and fMR data.
 *
 * This class performs the followeing steps:
 *    For i in VectorSize
 *      Extract VectorIndex Image
 *      Register with the fixed image
 *      Resample Vector Image
 *      Insert Resampled data into a new image
 *    End
 *
 */

template <class TInputImage, class TOutputImage>
class ITK_EXPORT VectorImageRegisterVersorScaleSkewFilter :
  public         ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef VectorImageRegisterVersorScaleSkewFilter Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>
    Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                                 InputImageType;
  typedef typename InputImageType::Pointer            InputImagePointer;
  typedef typename InputImageType::ConstPointer       InputImageConstPointer;
  typedef typename InputImageType::RegionType         InputImageRegionType;
  typedef typename InputImageType::SizeType           InputImageSizeType;
  typedef typename InputImageType::SpacingType        InputImageSpacingType;
  typedef typename InputImageType::PointType          InputImagePointType;
  typedef typename InputImageType::PixelType          InputImagePixelType;
  typedef typename InputImageType::DirectionType      InputImageDirectionType;
  typedef typename InputImageType::IndexType          InputImageIndexType;
  typedef typename InputImageType::InternalPixelType  InputImageValueType;
  typedef typename InputImagePixelType::ComponentType InputImageVectorComponentType;

#if 0
  /* The pixel type of FixedImageType should be based on the
   * InputImageValueType. However, the Mac compiler 4.0.1
   * does not like this. Currently this is hardcoded
   */
  typedef itk::Image<signed short, 3> FixedImageType;
#else
  //  typedef typename itk::Image<unsigned short,3>  FixedImageType;
  typedef typename itk::Image<InputImageVectorComponentType, 3> FixedImageType;
#endif
  typedef typename FixedImageType::Pointer   FixedImagePointer;
  typedef typename FixedImageType::PixelType FixedImagePixelType;

  typedef TOutputImage                                OutputImageType;
  typedef typename OutputImageType::Pointer           OutputImagePointer;
  typedef typename OutputImageType::ConstPointer      OutputImageConstPointer;
  typedef typename OutputImageType::RegionType        OutputImageRegionType;
  typedef typename OutputImageType::PixelType         OutputImagePixelType;
  typedef typename OutputImageType::IndexType         OutputImageIndexType;
  typedef typename OutputImageType::InternalPixelType OutputImageValueType;

  typedef itk::Image<OutputImageValueType, 3> CastImageType;

  /* Internally Used Typedefs */
  typedef itk::VectorIndexSelectionCastImageFilter<InputImageType, FixedImageType> VectorIndexFilterType;

  typedef itk::ScaleSkewVersor3DTransform<double>  TransformType;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::MattesMutualInformationImageToImageMetric<
      FixedImageType,
      FixedImageType>        MetricType;

  typedef itk::LinearInterpolateImageFunction<
      FixedImageType,
      double>         InterpolatorType;

  typedef itk::ImageRegistrationMethod<
      FixedImageType,
      FixedImageType>        RegistrationType;

  typedef itk::CenteredTransformInitializer<TransformType,
                                            FixedImageType,
                                            FixedImageType
                                            >  TransformInitializerType;
  typedef itk::ResampleImageFilter<
      FixedImageType,
      FixedImageType>    ResampleFilterType;

#if 0
  /* The output image type should be CastImageType for the
   * CastFilter. However, the Mac compiler 4.0.1
   * does not like this. Currently this is hardcoded to be
   * the FixedImageType to make it a scalar
   */
  typedef itk::CastImageFilter<
      FixedImageType,
      FixedImageType>    CastFilterType;
#else
  typedef itk::CastImageFilter<FixedImageType,
                               CastImageType>    CastFilterType;
#endif

  typedef typename VectorIndexFilterType::Pointer    VectorIndexFilterPointer;
  typedef typename TransformType::Pointer            TransformTypePointer;
  typedef typename TransformType::VersorType         VersorType;
  typedef typename VersorType::VectorType            VectorType;
  typedef typename MetricType::Pointer               MetricTypePointer;
  typedef typename OptimizerType::Pointer            OptimizerTypePointer;
  typedef typename OptimizerType::ParametersType     OptimizerParameterType;
  typedef typename OptimizerType::ScalesType         OptimizerScalesType;
  typedef typename InterpolatorType::Pointer         InterpolatorTypePointer;
  typedef typename RegistrationType::Pointer         RegistrationTypePointer;
  typedef typename TransformInitializerType::Pointer TransformInitializerTypePointer;
  typedef typename ResampleFilterType::Pointer       ResampleFilterTypePointer;
  typedef typename CastFilterType::Pointer           CastFilterTypePointer;

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

  /** The dimension of the input image must be 3 */
  itkConceptMacro( DimensionShouldBe3,
                   ( Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension), 3> ) );

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(VectorImageRegisterVersorScaleSkewFilter, ImageToImageFilter);

  /* SetInput and GetOutput Macros */
  itkSetObjectMacro(FixedImage, FixedImageType);
  itkGetConstObjectMacro(Output, OutputImageType);

  itkSetMacro(NumberOfSpatialSamples, int);
  itkSetMacro(NumberOfHistogramBins, int);
  itkSetMacro(NumberOfIterations, int);
  itkSetMacro(TranslationScale, float);
  itkSetMacro(ScalingScale, float);
  itkSetMacro(SkewScale, float);
  itkSetMacro(MaximumStepLength, float);
  itkSetMacro(MinimumStepLength, float);
  itkSetMacro(RelaxationFactor, float);
  itkSetMacro(OutputParameterFile, std::string);

  itkGetMacro(NumberOfSpatialSamples, int);
  itkGetMacro(NumberOfHistogramBins, int);
  itkGetMacro(NumberOfIterations, int);
  itkGetMacro(TranslationScale, float);
  itkGetMacro(ScalingScale, float);
  itkGetMacro(SkewScale, float);
  itkGetMacro(MaximumStepLength, float);
  itkGetMacro(MinimumStepLength, float);
  itkGetMacro(RelaxationFactor, float);
  itkGetMacro(OutputParameterFile, std::string);
protected:
  VectorImageRegisterVersorScaleSkewFilter();
  ~VectorImageRegisterVersorScaleSkewFilter()
  {
  }

  void GenerateData();

private:
  VectorImageRegisterVersorScaleSkewFilter(const Self &); // purposely not
                                                          // implemented
  void operator=(const Self &);                           // purposely not

  // implemented

  // Input and Output Image
  FixedImagePointer  m_FixedImage;
  OutputImagePointer m_Output;

  // Internally used Parameters
  float       m_TranslationScale;
  float       m_ScalingScale;
  float       m_SkewScale;
  float       m_MaximumStepLength;
  float       m_MinimumStepLength;
  float       m_RelaxationFactor;
  int         m_NumberOfSpatialSamples;
  int         m_NumberOfIterations;
  int         m_NumberOfHistogramBins;
  std::string m_OutputParameterFile;
};  // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorImageRegisterVersorScaleSkewFilter.hxx"
#endif

#endif
