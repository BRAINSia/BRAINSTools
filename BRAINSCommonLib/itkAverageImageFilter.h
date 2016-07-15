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
#ifndef __itkAverageImageFilter_h
#define __itkAverageImageFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"

namespace itk
{
/** \class AverageImageFilter
 *
 * \brief This filter performs pixelwise averaging among an arbitrary number
 * of input images.
 *
 * \par INPUTS
 * Input volumes must all contain the same size RequestedRegions. All input
 * images must have the same pixel type. All pixel types are supported that
 * supply operator=(), operator+=(), operator*=(double), and a cast to the
 * chosen output type (if not identical to input pixel type).
 *
 * \par OUTPUTS
 * The averaging filter produces a single output volume. Each output pixel
 * is assigned the average of the values from the corresponding input image
 * pixels.
 *
 * \par LIMITATIONS
 * For integer output pixel types, the result of the averaging is converted
 * to that type by a cast operation. There is currently no rounding
 * implemented.
 *
 * \author Torsten Rohlfing, SRI International, Neuroscience Program
 *
 * Funding for the implementation of this class was provided by NIAAA under
 * Grant No. AA05965, "CNS Deficits: Interaction of Age and Alcoholism",
 * PI: A. Pfefferbaum; and Grant No. AA13521, "INIA: Imaging Core",
 * PI: A. Pfefferbaum.
 */
template <typename TInputImage, typename TOutputImage = TInputImage>
class ITK_EXPORT AverageImageFilter :
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef AverageImageFilter Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(AverageImageFilter, ImageToImageFilter);

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef typename TInputImage::PixelType InputPixelType;

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Image typedef support */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;
  typedef typename InputImageType::ConstPointer InputImagePointer;
  typedef typename OutputImageType::Pointer     OutputImagePointer;

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;


protected:
  AverageImageFilter() {}
  virtual ~AverageImageFilter() {}

  void ThreadedGenerateData
  ( const OutputImageRegionType &outputRegionForThread, ThreadIdType threadId) ITK_OVERRIDE;

  void PrintSelf(std::ostream&, Indent) const ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(AverageImageFilter);
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAverageImageFilter.hxx"
#endif

#endif
