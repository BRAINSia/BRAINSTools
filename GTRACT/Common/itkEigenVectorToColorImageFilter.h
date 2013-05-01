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

#ifndef __itkEigenVectorToColorImageFilter_h
#define __itkEigenVectorToColorImageFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include "itkRGBAPixel.h"

#include <map>
#include <string>

namespace itk
{
/** \class TensorToAnisotropyImageFilter
 * \brief Calculates the Specified Anisotropy Index.
 *
 * The following Anisotropy Image are supported:
 *    Fractional Anistropy
 *    Relatibve Anisotropy
 *    Volume Ratio
 *    Radial Diffusivity
 *    Axial Diffusivity
 *    Coheernce Index
 *    Lattice Index
 *    Mean Diffusivity
 *
 */

enum ENUM_TENSOR_SHAPE_TYPE
  {
  TERTIARY_EIGENVECTOR = 0,
  SECONDARY_EIGENVECTOR = 1,
  PRIMARY_EIGENVECTOR = 2,
  TENSOR_SHAPE = 3
  };
typedef enum ENUM_TENSOR_SHAPE_TYPE TensorShapeType;

class EigenVectorToColorImageFilter : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef EigenVectorToColorImageFilter Self;
  typedef itk::Object                   Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Some convenient typedefs. */
  typedef itk::Vector<float, 6>         InputPixelType;
  typedef itk::Image<InputPixelType, 3> InputImageType;
  typedef InputImageType::Pointer       InputImagePointer;
  typedef InputImageType::ConstPointer  InputImageConstPointer;
  typedef InputImageType::RegionType    InputImageRegionType;
  typedef InputImageType::SizeType      InputImageSizeType;
  typedef InputImageType::SpacingType   InputImageSpacingType;
  typedef InputImageType::PointType     InputImagePointType;
  typedef InputImageType::PixelType     InputImagePixelType;
  typedef InputImageType::DirectionType InputImageDirectionType;

  typedef itk::RGBAPixel<unsigned char>  OutputPixelType;
  typedef itk::Image<OutputPixelType, 3> OutputImageType;
  typedef OutputImageType::Pointer       OutputImagePointer;
  typedef OutputImageType::RegionType    OutputImageRegionType;


/** Standard New method. */
  itkNewMacro(Self);

/** Runtime information support. */
  itkTypeMacro(EigenVectorToColorImageFilter, itk::Object);

  itkSetObjectMacro(Input,  InputImageType);
  itkGetConstObjectMacro(Output, OutputImageType);
  itkSetMacro(TensorShapeType, TensorShapeType);

  void Update();

protected:
  EigenVectorToColorImageFilter();
  ~EigenVectorToColorImageFilter()
  {
  }

private:
  EigenVectorToColorImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                // purposely not implemented

// Input and Output Image
  InputImagePointer  m_Input;
  OutputImagePointer m_Output;

  TensorShapeType m_TensorShapeType;
};  // end of class
} // end namespace itk

#endif
