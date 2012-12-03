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

#ifndef __itkGtractImageIO_h
#define __itkGtractImageIO_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkRGBAPixel.h"

#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaImageIO.h"
#include "gtractCommonWin32.h"

namespace itk
{
/** \class GtractImageIO
 * \brief Convience functions for image I/O. These were required
 * for building on Windows Visual studio because the object sizes
 * from the templated code.
 */
class GTRACT_COMMON_EXPORT GtractImageIO : public itk::Object
{
public:
  typedef GtractImageIO                 Self;
  typedef itk::Object                   Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(GtractImageIO, itk::Object);

  typedef itk::Image<signed short, 3>    Short3dImageType;
  typedef itk::Image<signed short, 4>    Short4dImageType;
  typedef itk::Image<float, 3>           Float3dImageType;
  typedef itk::RGBAPixel<unsigned char>  RGBPixelType;
  typedef itk::Image<RGBPixelType, 3>    Rgb3dImageType;
  typedef itk::Vector<float, 6>          TensorPixelType;
  typedef itk::Image<TensorPixelType, 3> TensorImageType;

  /*** Get/Set the Images for I/O Routines ***/
  itkSetObjectMacro(Short3dImage, Short3dImageType);
  itkSetObjectMacro(Short4dImage, Short4dImageType);
  itkSetObjectMacro(Float3dImage, Float3dImageType);
  itkSetObjectMacro(Rgb3dImage,   Rgb3dImageType);
  itkSetObjectMacro(TensorImage,  TensorImageType);

  itkGetConstObjectMacro(Short3dImage, Short3dImageType);
  itkGetConstObjectMacro(Short4dImage, Short4dImageType);
  itkGetConstObjectMacro(Float3dImage, Float3dImageType);
  itkGetConstObjectMacro(Rgb3dImage,   Rgb3dImageType);
  itkGetConstObjectMacro(TensorImage,  TensorImageType);

  /*** Additional API - Functions ***/
  void Load3dDICOMSeries();

  void Load4dDICOMSeries();

  void Load4dShortImage();

  void Load3dShortImage();

  void Load3dFloatImage();

  void Load3dRgbImage();

  void LoadTensorImage();

  void Save3dFloatImage();

  void Save3dShortImage();

  void Save4dShortImage();

  void Save3dRgbImage();

  void SaveTensorImage();

  void SetFileName(char *);

  void SetFileName(std::string);

  void SetDicomSeriesUID(char *);

  void SetDicomSeriesUID(std::string);

  void SetDicomDirectory(char *);

  void SetDicomDirectory(std::string);

protected:

  /** Constructor and Destructor */
  GtractImageIO();
  ~GtractImageIO()
  {
  }

private:
  GtractImageIO( const Self & );        // purposely not implemented
  void operator=( const Self & );       // purposely not implemented

  Short3dImageType::Pointer m_Short3dImage;
  Short4dImageType::Pointer m_Short4dImage;
  Float3dImageType::Pointer m_Float3dImage;
  Rgb3dImageType::Pointer   m_Rgb3dImage;
  TensorImageType::Pointer  m_TensorImage;

  /*** File / Directory / Dicom UID Strinmgs ***/
  std::string m_FileName;
  std::string m_DicomSeriesUID;
  std::string m_DicomDirectory;
}; // end of class
} // end namespace itk
#endif
