/*=========================================================================

  Program:   KWImage - Kitware Image IO Library
  Module:    $RCSfile: vtkKWImage.h,v $

  Copyright (c) Kitware, Inc., Insight Consortium.  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __vtkKWImage_h
#define __vtkKWImage_h

#include "itkImageBase.h"
#include "itkImageIOBase.h"

#include "vtkObject.h"

class vtkImageData;
class vtkImageImport;

/** \class KWImage
 *
 * \brief This class represents an image object.
 *
 * This class associates an internal ITK image and a VTK importer in such a way
 * that internally it can make available both image formats to ITK and VTK
 * classes.
 *
 */
class vtkKWImage : public vtkObject
{
public:
  static vtkKWImage * New();

  vtkTypeRevisionMacro(vtkKWImage, vtkObject);

  typedef itk::ImageBase<3>                 ImageBaseType;
  typedef ImageBaseType::Pointer            ImagePointer;
  typedef ImageBaseType::ConstPointer       ImageConstPointer;
  typedef itk::ImageIOBase::IOComponentType ITKScalarPixelType;

  // Set the untyped ITK image
  void SetITKImageBase( ImageBaseType * );

  // Return the pixel type using ITK enums.
  ITKScalarPixelType GetITKScalarPixelType() const;

  // Return the pixel type using VTK enums.
  int GetVTKScalarPixelType();

  vtkImageData * GetVTKImage();

  // Return the ITK image base. This is independent of the pixel type
  const ImageBaseType * GetITKImageBase() const;

protected:
  vtkKWImage();
  ~vtkKWImage();
private:
  vtkKWImage(const vtkKWImage &);      // Not implemented.
  void operator=(const vtkKWImage &);  // Not implemented.

  ImagePointer ItkImage;

  itk::ProcessObject::Pointer Exporter;

  vtkImageImport *Importer;
};

#endif
