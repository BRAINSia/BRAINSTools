/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVTKImageToImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2009/03/16 13:35:43 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkVTKImageToImageFilter_hxx
#define _itkVTKImageToImageFilter_hxx

#include "itkVTKImageToImageFilter.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TOutputImage>
VTKImageToImageFilter<TOutputImage>
::VTKImageToImageFilter()
{
  m_Exporter = vtkImageExport::New();

  m_Importer = ImporterFilterType::New();

  m_Importer->SetUpdateInformationCallback( m_Exporter->GetUpdateInformationCallback() );
  m_Importer->SetPipelineModifiedCallback( m_Exporter->GetPipelineModifiedCallback() );
  m_Importer->SetWholeExtentCallback( m_Exporter->GetWholeExtentCallback() );
  m_Importer->SetSpacingCallback( m_Exporter->GetSpacingCallback() );
  m_Importer->SetOriginCallback( m_Exporter->GetOriginCallback() );
  m_Importer->SetScalarTypeCallback( m_Exporter->GetScalarTypeCallback() );
  m_Importer->SetNumberOfComponentsCallback( m_Exporter->GetNumberOfComponentsCallback() );
  m_Importer->SetPropagateUpdateExtentCallback( m_Exporter->GetPropagateUpdateExtentCallback() );
  m_Importer->SetUpdateDataCallback( m_Exporter->GetUpdateDataCallback() );
  m_Importer->SetDataExtentCallback( m_Exporter->GetDataExtentCallback() );
  m_Importer->SetBufferPointerCallback( m_Exporter->GetBufferPointerCallback() );
  m_Importer->SetCallbackUserData( m_Exporter->GetCallbackUserData() );
}

/**
 * Destructor
 */
template <typename TOutputImage>
VTKImageToImageFilter<TOutputImage>
::~VTKImageToImageFilter()
{
  if( m_Exporter )
    {
    m_Exporter->Delete();
    m_Exporter = nullptr;
    }
}

/**
 * Set a vtkImageData as input
 */
template <typename TOutputImage>
void
VTKImageToImageFilter<TOutputImage>
::SetInputData( vtkImageData *inputImage )
{
  m_Exporter->SetInputData( inputImage );
}

/**
 * Get an itk::Image as output
 */
template <typename TOutputImage>
const typename VTKImageToImageFilter<TOutputImage>::OutputImageType
* VTKImageToImageFilter<TOutputImage>
::GetOutput() const
  {
  return m_Importer->GetOutput();
  }

/**
 * Get the exporter filter
 */
template <typename TOutputImage>
vtkImageExport *
VTKImageToImageFilter<TOutputImage>
::GetExporter() const
{
  return m_Exporter;
}

/**
 * Get the importer filter
 */
template <typename TOutputImage>
typename VTKImageToImageFilter<TOutputImage>::ImporterFilterType
* VTKImageToImageFilter<TOutputImage>
::GetImporter() const
  {
  return m_Importer;
  }

/**
 * Delegate the Update to the importer
 */
template <typename TOutputImage>
void
VTKImageToImageFilter<TOutputImage>
::Update()
{
  m_Importer->Update();
}

} // end namespace itk

#endif
