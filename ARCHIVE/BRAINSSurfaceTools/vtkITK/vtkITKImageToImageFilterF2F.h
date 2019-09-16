/*=========================================================================

  Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   vtkITK
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/

#ifndef __vtkITKImageToImageFilterF2F_h
#define __vtkITKImageToImageFilterF2F_h


#include "vtkITKImageToImageFilter.h"
#include "vtkImageAlgorithm.h"
#include "itkImageToImageFilter.h"
#include "itkVTKImageExport.h"
#include "itkVTKImageImport.h"
#include "vtkITKUtility.h"


class VTK_ITK_EXPORT vtkITKImageToImageFilterF2F : public vtkITKImageToImageFilter
{
public:
  vtkTypeMacro(vtkITKImageToImageFilterF2F, vtkITKImageToImageFilter);
  static vtkITKImageToImageFilterF2F *
  New()
  {
    return 0;
  };
  void
  PrintSelf(ostream & os, vtkIndent indent)
  {
    Superclass::PrintSelf(os, indent);
    os << m_Filter;
  };

protected:
  /// To/from ITK
  using InputImagePixelType = itk::Vector<float, 2>;
  using OutputImagePixelType = float;
  using InputImageType = itk::Image<InputImagePixelType, 3>;
  using OutputImageType = itk::Image<OutputImagePixelType, 3>;

  using ImageImportType = itk::VTKImageImport<InputImageType>;
  using ImageExportType = itk::VTKImageExport<OutputImageType>;
  ImageImportType::Pointer itkImporter;
  ImageExportType::Pointer itkExporter;

  using GenericFilterType = itk::ImageToImageFilter<InputImageType, OutputImageType>;
  GenericFilterType::Pointer m_Filter;

  vtkITKImageToImageFilterF2F(GenericFilterType * filter)
  {
    /// Need an import, export, and a ITK pipeline
    m_Filter = filter;
    this->itkImporter = ImageImportType::New();
    this->itkExporter = ImageExportType::New();
    ConnectPipelines(this->vtkExporter, this->itkImporter);
    ConnectPipelines(this->itkExporter, this->vtkImporter);
    this->LinkITKProgressToVTKProgress(m_Filter);

    /// Set up the filter pipeline
    m_Filter->SetInput(this->itkImporter->GetOutput());
    this->itkExporter->SetInput(m_Filter->GetOutput());
    this->vtkCast->SetOutputScalarTypeToFloat();
  };

  ~vtkITKImageToImageFilterF2F(){};

private:
  vtkITKImageToImageFilterF2F(const vtkITKImageToImageFilterF2F &); /// Not implemented.
  void
  operator=(const vtkITKImageToImageFilterF2F &); /// Not implemented.
};

#endif
