/*=========================================================================

  Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   vtkITK
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/

#ifndef __vtkITKImageToImageFilterFUL_h
#define __vtkITKImageToImageFilterFUL_h

#include "vtkITKImageToImageFilter.h"
#include "vtkImageAlgorithm.h"
#include "itkImageToImageFilter.h"
#include "itkVTKImageExport.h"
#include "itkVTKImageImport.h"
#include "vtkITKUtility.h"


class VTK_ITK_EXPORT vtkITKImageToImageFilterFUL : public vtkITKImageToImageFilter
{
public:
  vtkTypeMacro( vtkITKImageToImageFilterFUL, vtkITKImageToImageFilter );
  static vtkITKImageToImageFilterFUL *
  New()
  {
    return 0;
  };
  void
  PrintSelf( ostream & os, vtkIndent indent )
  {
    Superclass::PrintSelf( os, indent );
    os << m_Filter << std::endl;
  };

protected:
  /// To/from ITK
  using InputImageType = itk::Image< float, 3 >;
  using OutputImageType = itk::Image< unsigned long, 3 >;

  using ImageImportType = itk::VTKImageImport< InputImageType >;
  using ImageExportType = itk::VTKImageExport< OutputImageType >;
  ImageImportType::Pointer itkImporter;
  ImageExportType::Pointer itkExporter;

  using GenericFilterType = itk::ImageToImageFilter< InputImageType, OutputImageType >;
  GenericFilterType::Pointer m_Filter;

  vtkITKImageToImageFilterFUL( GenericFilterType * filter )
    : vtkITKImageToImageFilter()
  {
    /// Need an import, export, and a ITK pipeline
    m_Filter = filter;
    this->itkImporter = ImageImportType::New();
    this->itkExporter = ImageExportType::New();
    ConnectPipelines( this->vtkExporter, this->itkImporter );
    ConnectPipelines( this->itkExporter, this->vtkImporter );
    /// Set up the filter pipeline
    m_Filter->SetInput( this->itkImporter->GetOutput() );
    this->itkExporter->SetInput( m_Filter->GetOutput() );
    this->LinkITKProgressToVTKProgress( m_Filter );
    this->vtkCast->SetOutputScalarTypeToFloat();
  };

  ~vtkITKImageToImageFilterFUL(){};

private:
  vtkITKImageToImageFilterFUL( const vtkITKImageToImageFilterFUL & ); /// Not implemented.
  void
  operator=( const vtkITKImageToImageFilterFUL & ); /// Not implemented.
};

#endif
