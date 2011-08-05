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

#ifndef __itkFixDicomImageFilter_h
#define __itkFixDicomImageFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include "itkImageLinearConstIteratorWithIndex.h"

#include <map>
#include <string>

namespace itk
{
/** \class FixDicomImageFilter
 * \brief Fix the DICOM image. That is convert from 3D representation into 4D
 *
 * This class satisfies performs the following steps:
 *    For i in 4th Dimension
 *      ExtractVolume with Extract Image Filter
 *      Orient 3D extracted volume
 *    End
 *
 * It is build upon the ExtractImageFilter and the OrientImageFilter
 */

template <class TInputImage, class TOutputImage>
class FixDicomImageFilter : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef FixDicomImageFilter      Self;
  typedef itk::Object              Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                                            InputImageType;
  typedef TOutputImage                                           OutputImageType;
  typedef typename InputImageType::Pointer                       InputImagePointer;
  typedef typename InputImageType::ConstPointer                  InputImageConstPointer;
  typedef typename InputImageType::RegionType                    InputImageRegionType;
  typedef typename InputImageType::SizeType                      InputImageSizeType;
  typedef typename InputImageType::SpacingType                   InputImageSpacingType;
  typedef typename InputImageType::PointType                     InputImagePointType;
  typedef typename InputImageType::PixelType                     InputImagePixelType;
  typedef typename InputImageType::IndexType                     InputImageIndexType;
  typedef typename InputImageType::DirectionType                 InputImageDirectionType;
  typedef itk::ImageLinearConstIteratorWithIndex<InputImageType> IteratorType;

  typedef typename OutputImageType::Pointer       OutputImagePointer;
  typedef typename OutputImageType::ConstPointer  OutputImageConstPointer;
  typedef typename OutputImageType::RegionType    OutputImageRegionType;
  typedef typename OutputImageType::SizeType      OutputImageSizeType;
  typedef typename OutputImageType::SpacingType   OutputImageSpacingType;
  typedef typename OutputImageType::PointType     OutputImagePointType;
  typedef typename OutputImageType::PixelType     OutputImagePixelType;
  typedef typename OutputImageType::IndexType     OutputImageIndexType;
  typedef typename OutputImageType::DirectionType OutputImageDirectionType;

  /** ImageDimension constants * /
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** The dimensions of the input image must equal those of the
      output image. * /
  itkConceptMacro(SameDimension,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),itkGetStaticConstMacro(OutputImageDimension)>));

  /** The dimension of the input image must be 4. * /
  itkConceptMacro(DimensionShouldBe4,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),4>));
*/
  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(Orient4dImageFilter, itk::Object);

  /* SetInput and GetOutput Macros */
  itkSetObjectMacro(Input, InputImageType);
  itkGetObjectMacro(Output, OutputImageType);
  itkSetMacro(NumberOfImagesPerSlice, int);
  itkGetMacro(NumberOfImagesPerSlice, int);

  itkSetMacro(NumberBValues, int);
  itkGetMacro(NumberBValues, int);

  itkSetMacro(NumberDtiDirections, int);
  itkGetMacro(NumberDtiDirections, int);

  void SetDirectionsIncrementFirst();

  void SetBvaluesIncrementFirst();

  void Update();

protected:
  FixDicomImageFilter();
  ~FixDicomImageFilter()
  {
  }

private:
  FixDicomImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);      // purposely not implemented

  // Input and Output Image
  InputImagePointer  m_Input;
  OutputImagePointer m_Output;

  int m_NumberOfImagesPerSlice;
  int m_NumberBValues;
  int m_NumberDtiDirections;

  bool m_DirectionsIncrementFirst;
};  // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFixDicomImageFilter.txx"
#endif

#endif
