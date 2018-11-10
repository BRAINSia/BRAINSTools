/*=========================================================================

  Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   vtkITK
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/

#ifndef __vtkITKImageToImageFilterUSUL_h
#define __vtkITKImageToImageFilterUSUL_h

#include "vtkImageAlgorithm.h"
#include "vtkITKImageToImageFilter.h"
#include "itkImageToImageFilter.h"
#include "itkVTKImageExport.h"
#include "itkVTKImageImport.h"
#include "vtkITKUtility.h"


class VTK_ITK_EXPORT vtkITKImageToImageFilterUSUL : public vtkITKImageToImageFilter
{
public:
  vtkTypeMacro(vtkITKImageToImageFilterUSUL,vtkITKImageToImageFilter);
  static vtkITKImageToImageFilterUSUL* New() { return 0; };
  void PrintSelf(ostream& os, vtkIndent indent)
  {
    Superclass::PrintSelf ( os, indent );
    os << m_Filter;
  };

protected:

  /// To/from ITK
  using InputImageType = itk::Image<unsigned short, 3>;
  using OutputImageType = itk::Image<unsigned long, 3>;

  using ImageImportType = itk::VTKImageImport<InputImageType>;
  using ImageExportType = itk::VTKImageExport<OutputImageType>;
  ImageImportType::Pointer itkImporter;
  ImageExportType::Pointer itkExporter;

  using FilterType = itk::ImageToImageFilter<InputImageType,OutputImageType>;
  FilterType::Pointer m_Filter;

  vtkITKImageToImageFilterUSUL ( FilterType* filter )
  {
    /// Need an import, export, and a ITK pipeline
    m_Filter = filter;
    this->itkImporter = ImageImportType::New();
    this->itkExporter = ImageExportType::New();
    ConnectPipelines(this->vtkExporter, this->itkImporter);
    ConnectPipelines(this->itkExporter, this->vtkImporter);
    this->LinkITKProgressToVTKProgress ( m_Filter );
    /// Set up the filter pipeline
    m_Filter->SetInput ( this->itkImporter->GetOutput() );
    this->itkExporter->SetInput ( m_Filter->GetOutput() );
    this->vtkCast->SetOutputScalarTypeToUnsignedShort();
  };

  ~vtkITKImageToImageFilterUSUL()
  {
  };

private:
  vtkITKImageToImageFilterUSUL(const vtkITKImageToImageFilterUSUL&);  /// Not implemented.
  void operator=(const vtkITKImageToImageFilterUSUL&);  /// Not implemented.
};

#endif
