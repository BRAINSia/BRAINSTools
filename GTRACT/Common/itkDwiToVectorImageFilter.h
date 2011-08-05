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

#ifndef __itkDwiToVectorImageFilter_h
#define __itkDwiToVectorImageFilter_h

#include <itkObject.h>
#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkIOCommon.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkVariableLengthVector.h>
#include <itkOrientImageFilter.h>

#include <map>
#include <string>

namespace itk
{
/** \class FixDicomImageFilter
 * \brief Fix the DICOM image. That is convert from 3D representation into 4D
 *
 * This class performs the following steps:
 *    For i in 4th Dimension
 *      ExtractVolume with Extract Image Filter
 *      Rotate Gradients if specified by the user
 *    End
 *
 * It is build upon the ExtractImageFilter
 */

template <class TInputImage, class TOutputImage>
class DwiToVectorImageFilter : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef DwiToVectorImageFilter   Self;
  typedef itk::Object              Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::RegionType      InputImageRegionType;
  typedef typename InputImageType::SizeType        InputImageSizeType;
  typedef typename InputImageType::SpacingType     InputImageSpacingType;
  typedef typename InputImageType::PointType       InputImagePointType;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef typename InputImageType::IndexType       InputImageIndexType;
  typedef typename InputImageType::DirectionType   InputImageDirectionType;
  typedef ImageRegionConstIterator<InputImageType> ConstIteratorType;
  typedef ImageRegionIterator<InputImageType>      IteratorType;

  typedef TOutputImage                            OutputImageType;
  typedef typename OutputImageType::Pointer       OutputImagePointer;
  typedef typename OutputImageType::ConstPointer  OutputImageConstPointer;
  typedef typename OutputImageType::RegionType    OutputImageRegionType;
  typedef typename OutputImageType::SizeType      OutputImageSizeType;
  typedef typename OutputImageType::SpacingType   OutputImageSpacingType;
  typedef typename OutputImageType::PointType     OutputImagePointType;
  typedef typename OutputImageType::PixelType     OutputImagePixelType;
  typedef typename OutputImageType::IndexType     OutputImageIndexType;
  typedef typename OutputImageType::DirectionType OutputImageDirectionType;
  typedef typename OutputImageType::ValueType     OutputImageValueType;
  typedef ImageRegionIterator<OutputImageType>    VectorIteratorType;

  typedef SpatialOrientation::ValidCoordinateOrientationFlags CoordinateType;
  typedef vnl_matrix<float>                                   MatrixType;
  typedef vnl_vector<float>                                   VectorType;

  /*const unsigned int ImageDimension = TInputImage::ImageDimension;*/

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** The dimensions of the input image must equal those of the
      output image. * /
  itkConceptMacro(SameDimension,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),itkGetStaticConstMacro(OutputImageDimension)>));

  / ** The dimension of the input image must be 4. * /
  itkConceptMacro(DimensionShouldBe4,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),4>));
*/
  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(DwiToVectorImageFilter, itk::Object);

  /** SetInput and GetOutput Macros */
  itkSetObjectMacro(Input, InputImageType);
  itkGetObjectMacro(Output, OutputImageType);

  /** Image Sorting Function Macros */
  itkSetMacro(NumberDtiDirections, int);
  itkGetMacro(NumberDtiDirections, int);
  itkSetMacro(BValue, float);
  itkGetMacro(BValue, float);
  itkGetMacro(DiffusionDirections, MatrixType);
  itkSetMacro(FlipZ, bool);
  itkSetMacro(FlipY, bool);
  itkSetMacro(FlipX, bool);

  /** Image Orientation Function Macros */
  itkSetMacro(RotateGradients, bool);

  /** Image Order has DWI gradient directions incrementing faster than B values
    */
  void SetDirectionsIncrementFirst();

  /** Image Order has DWI B values incrementing faster than Gradient Directions
    */
  void SetBvaluesIncrementFirst();

  /** Set the Diffusion Directions */
  void SetDiffusionDirections( MatrixType );

  void Update();

protected:
  DwiToVectorImageFilter();
  ~DwiToVectorImageFilter()
  {
  }

  /** Set Image NRRD Meta Data - Used for Image I/O */
  void SetMetaDataHeader();

private:
  DwiToVectorImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);         // purposely not implemented

  // Input and Output Image
  InputImagePointer  m_Input;
  OutputImagePointer m_Output;

  // Image Sorting Information
  int m_NumberDtiDirections;

  // Diffusion Encoding Information
  MatrixType m_DiffusionDirections;
  double     m_BValue;

  // Orientation Information
  bool m_RotateGradients;

  // Flip Orientation Information
  bool m_FlipZ;
  bool m_FlipY;
  bool m_FlipX;
};  // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDwiToVectorImageFilter.txx"
#endif

#endif
