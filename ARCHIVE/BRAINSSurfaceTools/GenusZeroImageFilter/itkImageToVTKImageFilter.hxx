/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToVTKImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2009/03/16 13:35:43 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkImageToVTKImageFilter_hxx
#define _itkImageToVTKImageFilter_hxx

#include "itkImageToVTKImageFilter.h"

namespace itk
{
/**
 * Constructor
 */
template < typename TInputImage >
ImageToVTKImageFilter< TInputImage >::ImageToVTKImageFilter()
{
  m_Importer = vtkImageImport::New();

  m_Exporter = ExporterFilterType::New();

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
template < typename TInputImage >
ImageToVTKImageFilter< TInputImage >::~ImageToVTKImageFilter()
{
  if ( m_Importer )
  {
    m_Importer->Delete();
    m_Importer = nullptr;
  }
}

/**
 * Set an itk::Image as input
 */
template < typename TInputImage >
void
ImageToVTKImageFilter< TInputImage >::SetInputData( const InputImageType * inputImage )
{
  m_Exporter->SetInput( inputImage );
}

/**
 * Get a vtkImage as output
 */
template < typename TInputImage >
vtkImageData *
ImageToVTKImageFilter< TInputImage >::GetOutput() const
{
  return m_Importer->GetOutput();
}

/**
 * Get the importer filter
 */
template < typename TInputImage >
vtkImageImport *
ImageToVTKImageFilter< TInputImage >::GetImporter() const
{
  return m_Importer;
}

/**
 * Get the exporter filter
 */
template < typename TInputImage >
typename ImageToVTKImageFilter< TInputImage >::ExporterFilterType *
ImageToVTKImageFilter< TInputImage >::GetExporter() const
{
  return m_Exporter.GetPointer();
}

/**
 * Delegate the Update to the importer
 */
template < typename TInputImage >
void
ImageToVTKImageFilter< TInputImage >::Update()
{
  m_Importer->Update();
}

} // end namespace itk

#endif
