/*=========================================================================

  Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   vtkITK
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/

#ifndef __vtkITKImageToImageFilterUSUS_h
#define __vtkITKImageToImageFilterUSUS_h


#include "vtkITKImageToImageFilter.h"
#include "vtkImageAlgorithm.h"
#include "itkImageToImageFilter.h"
#include "itkVTKImageExport.h"
#include "itkVTKImageImport.h"
#include "vtkITKUtility.h"


class VTK_ITK_EXPORT vtkITKImageToImageFilterUSUS : public vtkITKImageToImageFilter
{
public:
  vtkTypeMacro( vtkITKImageToImageFilterUSUS, vtkITKImageToImageFilter );
  static vtkITKImageToImageFilterUSUS *
  New()
  {
    return 0;
  };
  void
  PrintSelf( ostream & os, vtkIndent indent )
  {
    Superclass::PrintSelf( os, indent );
    os << m_Filter;
  };

protected:
  /// To/from ITK
  using InputImagePixelType = unsigned short;
  using OutputImagePixelType = unsigned short;
  using InputImageType = itk::Image< InputImagePixelType, 3 >;
  using OutputImageType = itk::Image< OutputImagePixelType, 3 >;

  using ImageImportType = itk::VTKImageImport< InputImageType >;
  using ImageExportType = itk::VTKImageExport< OutputImageType >;
  ImageImportType::Pointer itkImporter;
  ImageExportType::Pointer itkExporter;

  using GenericFilterType = itk::ImageToImageFilter< InputImageType, OutputImageType >;
  GenericFilterType::Pointer m_Filter;

  vtkITKImageToImageFilterUSUS( GenericFilterType * filter )
  {
    /// Need an import, export, and a ITK pipeline
    m_Filter = filter;
    this->itkImporter = ImageImportType::New();
    this->itkExporter = ImageExportType::New();
    ConnectPipelines( this->vtkExporter, this->itkImporter );
    ConnectPipelines( this->itkExporter, this->vtkImporter );
    this->LinkITKProgressToVTKProgress( m_Filter );

    /// Set up the filter pipeline
    m_Filter->SetInput( this->itkImporter->GetOutput() );
    this->itkExporter->SetInput( m_Filter->GetOutput() );
    this->vtkCast->SetOutputScalarTypeToUnsignedShort();
  };

  ~vtkITKImageToImageFilterUSUS(){};

private:
  vtkITKImageToImageFilterUSUS( const vtkITKImageToImageFilterUSUS & ); /// Not implemented.
  void
  operator=( const vtkITKImageToImageFilterUSUS & ); /// Not implemented.
};

#endif
