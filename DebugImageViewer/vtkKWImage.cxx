/*=========================================================================

  Program:   KWImage - Kitware Image IO Library
  Module:    $RCSfile: vtkKWImage.cxx,v $

  Copyright (c) Kitware, Inc., Insight Consortium.  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkKWImage.h"

#include "itkImage.h"
#include "itkVTKImageExport.h"

#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkImageImport.h"

/** \class PipelineCreator
 *  This helper class will take care of instantiating the appropriate
 *  ITK Export class corresponding to the actual pixel type of the
 *  input image. */
template <class TPixel>
class PipelineCreator
{
public:

  typedef itk::ImageBase<3>           ImageBaseType;
  typedef ImageBaseType::Pointer      ImageBasePointer;
  typedef itk::ProcessObject          ExporterBaseType;
  typedef itk::ProcessObject::Pointer ExporterBasePointer;
  typedef itk::Image<TPixel, 3>       ImageType;

  static void CreateExporter( ImageBasePointer    & imageBase,
                              ExporterBasePointer & exporter,
                              vtkImageImport      *importer  )
  {
    ImageType *image
      = dynamic_cast<ImageType *>( imageBase.GetPointer() );

    if( image )
      {
      typedef itk::VTKImageExport<ImageType>     ExportFilterType;
      typedef typename ExportFilterType::Pointer ExportFilterPointer;
      ExportFilterPointer itkExporter = ExportFilterType::New();
      itkExporter->SetInput( image );

      exporter = itkExporter;

      importer->SetUpdateInformationCallback(
        itkExporter->GetUpdateInformationCallback() );
      importer->SetPipelineModifiedCallback(
        itkExporter->GetPipelineModifiedCallback() );
      importer->SetWholeExtentCallback(
        itkExporter->GetWholeExtentCallback() );
      importer->SetSpacingCallback(
        itkExporter->GetSpacingCallback() );
      importer->SetOriginCallback(
        itkExporter->GetOriginCallback() );
      importer->SetScalarTypeCallback(
        itkExporter->GetScalarTypeCallback() );
      importer->SetNumberOfComponentsCallback(
        itkExporter->GetNumberOfComponentsCallback() );
      importer->SetPropagateUpdateExtentCallback(
        itkExporter->GetPropagateUpdateExtentCallback() );
      importer->SetUpdateDataCallback(
        itkExporter->GetUpdateDataCallback() );
      importer->SetDataExtentCallback(
        itkExporter->GetDataExtentCallback() );
      importer->SetBufferPointerCallback(
        itkExporter->GetBufferPointerCallback() );
      importer->SetCallbackUserData(
        itkExporter->GetCallbackUserData() );
      }
  }
};

/** This helper macro will instantiate the pipeline creator for a particular
 * pixel type */
#define CreatePipelineMacro( PixelType )                \
  PipelineCreator<PixelType>::CreateExporter(          \
    this->ItkImage, this->Exporter, this->Importer );

// ----------------------------------------------------------------------------
vtkStandardNewMacro( vtkKWImage );
vtkCxxRevisionMacro(vtkKWImage, "$Revision: 1.1 $");

// ----------------------------------------------------------------------------
vtkKWImage::vtkKWImage()
{
  this->Importer = vtkImageImport::New();
}

// ----------------------------------------------------------------------------
vtkKWImage::~vtkKWImage()
{
  if( this->Importer )
    {
    this->Importer->Delete();
    }
}

// ----------------------------------------------------------------------------
void vtkKWImage::SetITKImageBase( ImageBaseType *image )
{
  if( this->ItkImage.GetPointer() == image )
    {
    return;
    }

  this->ItkImage = image;
  this->Modified();

  CreatePipelineMacro( unsigned char );
  CreatePipelineMacro( char );
  CreatePipelineMacro( unsigned short );
  CreatePipelineMacro( short );
  CreatePipelineMacro( unsigned int );
  CreatePipelineMacro( int );
  CreatePipelineMacro( unsigned long );
  CreatePipelineMacro( long );
  CreatePipelineMacro( float );
  CreatePipelineMacro( double );

  this->Importer->Update();
}

// ----------------------------------------------------------------------------
vtkImageData * vtkKWImage::GetVTKImage()
{
  return this->Importer->GetOutput();
}

// ----------------------------------------------------------------------------
const vtkKWImage::ImageBaseType * vtkKWImage::GetITKImageBase() const
{
  return this->ItkImage;
}

// ----------------------------------------------------------------------------
vtkKWImage::ITKScalarPixelType vtkKWImage::GetITKScalarPixelType() const
{
  ITKScalarPixelType pixelType = itk::ImageIOBase::UCHAR;

  ImageBaseType *itkImageBase = this->ItkImage.GetPointer();

  if( dynamic_cast<itk::Image<unsigned char, 3> *>( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::UCHAR;
    }
  else if( dynamic_cast<itk::Image<char, 3> *>( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::CHAR;
    }
  else if( dynamic_cast<itk::Image<short, 3> *>( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::SHORT;
    }
  else if( dynamic_cast<itk::Image<unsigned short, 3> *>( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::USHORT;
    }
  else if( dynamic_cast<itk::Image<int, 3> *>( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::INT;
    }
  else if( dynamic_cast<itk::Image<unsigned int, 3> *>( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::UINT;
    }
  else if( dynamic_cast<itk::Image<long, 3> *>( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::LONG;
    }
  else if( dynamic_cast<itk::Image<unsigned long, 3> *>( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::ULONG;
    }
  else if( dynamic_cast<itk::Image<float, 3> *>( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::FLOAT;
    }
  else if( dynamic_cast<itk::Image<double, 3> *>( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::DOUBLE;
    }

  return pixelType;
}

// ----------------------------------------------------------------------------
int vtkKWImage::GetVTKScalarPixelType()
{
  return this->GetVTKImage()->GetScalarType();
}
