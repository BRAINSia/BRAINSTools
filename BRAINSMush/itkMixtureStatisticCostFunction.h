/*=========================================================================
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile: itkMixtureStatisticsCostFunction.h,v $
 *  Language:  C++
 *  Date:      $Date: 2007-03-29 19:37:00 $
 *  Version:   $Revision: 1.14 $
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef __itkMixtureStatisticsCostFunction_h
#define __itkMixtureStatisticsCostFunction_h

#include "itkImageBase.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkMultipleValuedCostFunction.h"
#include "itkExceptionObject.h"

namespace itk
{
template <class TFirstImage, class TSecondImage>
class MixtureStatisticCostFunction :
  public MultipleValuedCostFunction
{
public:

  /** Standard typedefs. */
  typedef MixtureStatisticCostFunction Self;
  typedef MultipleValuedCostFunction   Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(MixtureStatisticCostFunction, MultipleValuedCostFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /**  Type of the First Image. */
  typedef TFirstImage                           FirstImageType;
  typedef typename TFirstImage::PixelType       FirstImagePixelType;
  typedef typename FirstImageType::ConstPointer FirstImageConstPointer;

  /**  Type of the Second Image. */
  typedef TSecondImage                           SecondImageType;
  typedef typename TFirstImage::PixelType        SecondImagePixelType;
  typedef typename SecondImageType::ConstPointer SecondImageConstPointer;
  typedef typename SecondImageType::RegionType   SecondImageRegionType;

  /** Constants for the image dimensions */
  itkStaticConstMacro(FirstImageDimension,
                      unsigned int,
                      TFirstImage::ImageDimension);
  itkStaticConstMacro(SecondImageDimension,
                      unsigned int,
                      TSecondImage::ImageDimension);

  /** Array Typedefs. */
  typedef Superclass::ParametersType ParametersType;
  typedef Superclass::MeasureType    MeasureType;
  typedef Superclass::DerivativeType DerivativeType;

  /**  Type for the mask of the first image. Only pixels that are "inside"
    * this mask will be considered for the computation of the metric */
  typedef itk::Image<signed short, 3>      ImageMaskType;
  typedef typename  ImageMaskType::Pointer ImageMaskPointer;

  /** Connect the First Image.  */
  itkSetConstObjectMacro(FirstImage, FirstImageType);

  /** Get the First Image. */
  itkGetConstObjectMacro(FirstImage, FirstImageType);

  /** Connect the Second Image.  */
  itkSetConstObjectMacro(SecondImage, SecondImageType);

  /** Get the Second Image. */
  itkGetConstObjectMacro(SecondImage, SecondImageType);

  /** Set/Get the first image mask. */
  itkSetObjectMacro(ImageMask, ImageMaskType);
  itkGetConstObjectMacro(ImageMask, ImageMaskType);

  itkGetMacro(DesiredMean, double);
  itkSetMacro(DesiredMean, double);
  itkGetMacro(DesiredVariance, double);
  itkSetMacro(DesiredVariance, double);

  itkGetMacro(NumberOfMaskVoxels, double);
  itkGetMacro(SumOfFirstMaskVoxels, double);
  itkGetMacro(SumOfSecondMaskVoxels, double);
  itkGetMacro(SumSquaresOfFirstMaskVoxels, double);
  itkGetMacro(SumSquaresOfSecondMaskVoxels, double);
  itkGetMacro(SumOfFirstTimesSecondMaskVoxels, double);

  /** The dimensions of parameter space. */
  enum { SpaceDimension = 2 };

  /** Not necessary for this optimizer. */
  void GetDerivative( const ParametersType & itkNotUsed(parameters),
                      DerivativeType & itkNotUsed(derivative) ) const
  {
  }

  /** Return the values evaluated for the given parameters. */
  MeasureType GetValue(const ParametersType & parameters) const;

  /** Return a pointer of values evaluated for the given parameters. */
  MeasureType * GetValue(ParametersType & parameters);

  /** Get the SpaceDimension. */
  unsigned int GetNumberOfParameters() const;

  /** Get the number Range Dimension. */
  unsigned int GetNumberOfValues() const;

  /** Initialize */
  void Initialize(short label);

protected:
  MixtureStatisticCostFunction();
  virtual ~MixtureStatisticCostFunction();

  void PrintSelf(std::ostream & os, Indent indent) const;

  FirstImageConstPointer  m_FirstImage;
  SecondImageConstPointer m_SecondImage;

  mutable ImageMaskPointer m_ImageMask;
private:
  MixtureStatisticCostFunction(const Self &);   // purposely not implemented
  void operator=(const Self &);                 // purposely not implemented

  double m_DesiredMean;
  double m_DesiredVariance;

  double m_NumberOfMaskVoxels;
  double m_SumOfFirstMaskVoxels;
  double m_SumOfSecondMaskVoxels;
  double m_SumSquaresOfFirstMaskVoxels;
  double m_SumSquaresOfSecondMaskVoxels;
  double m_SumOfFirstTimesSecondMaskVoxels;

  /** Different arrays. */
  mutable MeasureType    m_Measure;
  mutable MeasureType *  m_MeasurePointer;
  mutable ParametersType m_Parameters;
};
}   // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMixtureStatisticCostFunction.hxx"
#endif

#endif
