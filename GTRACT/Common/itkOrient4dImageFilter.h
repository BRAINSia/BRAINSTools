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

#ifndef __itkOrient4dImageFilter_h
#define __itkOrient4dImageFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkSpatialOrientation.h"
#include "itkOrientImageFilter.h"
#include "itkExtractImageFilter.h"

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

template <class TInputImage, class TOutputImage>
class Orient4dImageFilter : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef Orient4dImageFilter      Self;
  typedef itk::Object              Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                            InputImageType;
  typedef TOutputImage                           OutputImageType;
  typedef typename InputImageType::Pointer       InputImagePointer;
  typedef typename InputImageType::ConstPointer  InputImageConstPointer;
  typedef typename InputImageType::RegionType    InputImageRegionType;
  typedef typename InputImageType::SizeType      InputImageSizeType;
  typedef typename InputImageType::SpacingType   InputImageSpacingType;
  typedef typename InputImageType::PointType     InputImagePointType;
  typedef typename InputImageType::PixelType     InputImagePixelType;
  typedef typename InputImageType::DirectionType InputImageDirectionType;

  typedef typename OutputImageType::Pointer      OutputImagePointer;
  typedef typename OutputImageType::ConstPointer OutputImageConstPointer;
  typedef typename OutputImageType::RegionType   OutputImageRegionType;
  typedef typename OutputImageType::PixelType    OutputImagePixelType;
  typedef typename OutputImageType::IndexType    OutputImageIndexType;

  typedef itk::Image<InputImagePixelType, 3>       ExtractImageType;
  typedef typename ExtractImageType::Pointer       ExtractImagePointer;
  typedef typename ExtractImageType::DirectionType ExtractImageDirectionType;
  typedef itk::Image<OutputImagePixelType, 3>      OrientImageType;
  typedef typename OrientImageType::Pointer        OrientImagePointer;
  typedef typename OrientImageType::IndexType      OrientImageIndexType;
  typedef typename OrientImageType::RegionType     OrientImageRegionType;
  typedef typename OrientImageType::SpacingType    OrientImageSpacingType;
  typedef typename OrientImageType::PointType      OrientImagePointType;
  typedef typename OrientImageType::SizeType       OrientImageSizeType;

  typedef SpatialOrientation::ValidCoordinateOrientationFlags
    CoordinateOrientationCode;
  /** Axes permuter type. */
  typedef PermuteAxesImageFilter<TInputImage>          PermuterType;
  typedef typename PermuterType::PermuteOrderArrayType PermuteOrderArrayType;

  /** Axes flipper type. */

  typedef ExtractImageFilter<TInputImage, ExtractImageType>    ExtractFilterType;
  typedef typename ExtractFilterType::Pointer                  ExtractFilterPointerType;
  typedef OrientImageFilter<ExtractImageType, OrientImageType> OrientFilterType;
  typedef typename OrientFilterType::Pointer                   OrientFilterPointerType;
  typedef FlipImageFilter<OrientImageType>                     FlipFilterType;
  typedef typename FlipFilterType::Pointer                     FlipFilterPointerType;
  typedef typename FlipFilterType::FlipAxesArrayType           FlipFilterAxesType;

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

  itkSetMacro(FlipXaxis, int);
  itkSetMacro(FlipYaxis, int);
  itkSetMacro(FlipZaxis, int);

  /**********************************************************************
   * Methods for Setting the Resulting Orientation to a particular plane
   ***********************************************************************/
  void SetDesiredCoordinateOrientationToAxial()
  {
    m_OrientImageFilter->SetDesiredCoordinateOrientationToAxial();
  }

  void SetDesiredCoordinateOrientationToCoronal()
  {
    m_OrientImageFilter->SetDesiredCoordinateOrientationToCoronal();
  }

  void SetDesiredCoordinateOrientationToSagittal()
  {
    m_OrientImageFilter->SetDesiredCoordinateOrientationToSagittal();
  }

  /**********************************************************************
   * Turn On/Off use of direction vcl_cosines in the image
   ***********************************************************************/
  void UseImageDirectionOn()
  {
    m_OrientImageFilter->UseImageDirectionOn();
  }

  void UseImageDirectionOff()
  {
    m_OrientImageFilter->UseImageDirectionOff();
  }

  /**********************************************************************
   * Get/Set the Orientation values in the image
   ***********************************************************************/

  void SetGivenCoordinateOrientation(CoordinateOrientationCode newCode)
  {
    m_OrientImageFilter->SetGivenCoordinateOrientation(newCode);
  }

  void SetDesiredCoordinateOrientation(CoordinateOrientationCode newCode)
  {
    m_OrientImageFilter->SetDesiredCoordinateOrientation(newCode);
  }

  CoordinateOrientationCode GetGivenCoordinateOrientation()
  {
    return m_OrientImageFilter->GetGivenCoordinateOrientation();
  }

  CoordinateOrientationCode GetDesiredCoordinateOrientation()
  {
    return m_OrientImageFilter->GetDesiredCoordinateOrientation();
  }

  void Update();

protected:
  Orient4dImageFilter();
  ~Orient4dImageFilter()
  {
  }

private:
  Orient4dImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);      // purposely not implemented

  // Input and Output Image
  InputImagePointer  m_Input;
  OutputImagePointer m_Output;

  // Filters used internally by the program
  ExtractFilterPointerType m_ExtractImageFilter;
  OrientFilterPointerType  m_OrientImageFilter;

  // Optional Flip - Used to fix problems with the Direction Cosines
  int m_FlipXaxis;
  int m_FlipYaxis;
  int m_FlipZaxis;
};  // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOrient4dImageFilter.hxx"
#endif

#endif
