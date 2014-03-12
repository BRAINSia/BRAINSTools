/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
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

#ifndef __itkDwiToVectorImageFilter_hxx
#define __itkDwiToVectorImageFilter_hxx

#include <iostream>
#include <vector>

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkOrientImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkIOCommon.h>
#include <itkMetaDataObject.h>
#include <itkProgressAccumulator.h>
#include <itkFlipImageFilter.h>

#include "itkDwiToVectorImageFilter.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
DwiToVectorImageFilter<TInputImage, TOutputImage>
::DwiToVectorImageFilter()
{
  m_NumberDtiDirections = 1;
  m_RotateGradients = false;
  m_FlipZ = false;
  m_FlipY = false;
  m_FlipX = false;
}

template <class TInputImage, class TOutputImage>
void
DwiToVectorImageFilter<TInputImage, TOutputImage>
::SetDiffusionDirections( MatrixType diffusionDirections )
{
  m_DiffusionDirections = diffusionDirections;
}

template <class TInputImage, class TOutputImage>
void
DwiToVectorImageFilter<TInputImage, TOutputImage>
::Update()
{
  InputImageRegionType region3D  = m_Input->GetLargestPossibleRegion();
  InputImageIndexType  index3D  = region3D.GetIndex();
  InputImageSizeType   size3D    = region3D.GetSize();

  const int numberOfSlices = size3D[2] / m_NumberDtiDirections;

  size3D[2] = numberOfSlices;

  region3D.SetIndex(index3D);
  region3D.SetSize(size3D);

  typedef itk::Image<InputImagePixelType, InputImageDimension> TempImageType;
  typedef typename TempImageType::Pointer                      TempImagePointer;
  TempImagePointer tmpImage =  TempImageType::New();
  tmpImage->SetRegions( region3D );
  tmpImage->SetSpacing( m_Input->GetSpacing() );
  tmpImage->SetOrigin( m_Input->GetOrigin() );
  tmpImage->SetDirection( m_Input->GetDirection() );
  tmpImage->Allocate();

  typedef itk::Image<OutputImageValueType, OutputImageDimension> VectorIndexImageType;
  typedef typename VectorIndexImageType::Pointer                 VectorImageIndexPointer;
  VectorImageIndexPointer vectorIndexImage;

  IteratorType tempIt( tmpImage, tmpImage->GetRequestedRegion() );

  ConstIteratorType inputIt( m_Input, m_Input->GetRequestedRegion() );

  OutputImagePixelType vectorImagePixel( m_NumberDtiDirections );

  /* Should Progress Information be added ???? */
  int tmpIndex = 0;
  inputIt.GoToBegin();

  while( !inputIt.IsAtEnd() )
    {
    /* Replace with Extract Image Region Filter */
    for( tempIt.GoToBegin(); !tempIt.IsAtEnd(); ++tempIt, ++inputIt )
      {
      tempIt.Set( inputIt.Value() );
      }

    /* Flip Image if Specified */
    typedef FlipImageFilter<InputImageType> FlipImageFilterType;
    typename FlipImageFilterType::FlipAxesArrayType flipArray;
    flipArray.Fill(0);
    if( m_FlipZ )
      {
      flipArray[2] = 1;
      }
    if( m_FlipY )
      {
      flipArray[1] = 1;
      }
    if( m_FlipX )
      {
      flipArray[0] = 1;
      }

    typename FlipImageFilterType::Pointer flipImageFilter =  FlipImageFilterType::New();
    flipImageFilter->SetInput( tmpImage );
    flipImageFilter->SetFlipAxes( flipArray );
    flipImageFilter->Update();

    typedef CastImageFilter<InputImageType, VectorIndexImageType> CastFilterType;
    typedef typename CastFilterType::Pointer                      CastFilterPointer;
    CastFilterPointer CastImageFilter = CastFilterType::New();
    CastImageFilter->SetInput( flipImageFilter->GetOutput() );
    CastImageFilter->Update();
    vectorIndexImage = CastImageFilter->GetOutput();
    vectorIndexImage->SetDirection( m_Input->GetDirection() );

    if( tmpIndex == 0 )
      {
      std::cout << "Allocate Image: " << vectorIndexImage->GetLargestPossibleRegion() << std::endl;
      m_Output = OutputImageType::New();
      m_Output->SetRegions( vectorIndexImage->GetLargestPossibleRegion() );
      m_Output->SetSpacing( vectorIndexImage->GetSpacing() );
      m_Output->SetOrigin( vectorIndexImage->GetOrigin() );
      m_Output->SetDirection( vectorIndexImage->GetDirection() );
      m_Output->SetVectorLength( m_NumberDtiDirections );
      m_Output->Allocate();
      }

    typedef ImageRegionConstIterator<VectorIndexImageType> ConstVectorIndexIteratorType;
    ConstVectorIndexIteratorType vectorIndexIt( vectorIndexImage, vectorIndexImage->GetRequestedRegion() );

    VectorIteratorType outputIt( m_Output, m_Output->GetRequestedRegion() );
    for( outputIt.GoToBegin(), vectorIndexIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt, ++vectorIndexIt )
      {
      vectorImagePixel = outputIt.Get();
      vectorImagePixel[tmpIndex] = vectorIndexIt.Value();
      outputIt.Set( vectorImagePixel );
      }
    tmpIndex++;
    }

  SetMetaDataHeader();
}

template <class TInputImage, class TOutputImage>
void
DwiToVectorImageFilter<TInputImage, TOutputImage>
::SetMetaDataHeader()
{
  itk::MetaDataDictionary meta;

  std::string NrrdValue;

  NrrdValue = "cell";
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_centerings[0]", NrrdValue);
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_centerings[1]", NrrdValue);
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_centerings[2]", NrrdValue);
  NrrdValue = "none";
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_centerings[3]", NrrdValue);

  NrrdValue = "space";
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_kinds[0]", NrrdValue);
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_kinds[1]", NrrdValue);
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_kinds[2]", NrrdValue);
  NrrdValue = "list";
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_kinds[3]", NrrdValue);

  NrrdValue = "mm";
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_space units[0]", NrrdValue);
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_space units[1]", NrrdValue);
  itk::EncapsulateMetaData<std::string>(meta, "NRRD_space units[2]", NrrdValue);

  double spacing = ( m_Output->GetSpacing() )[2];
  itk::EncapsulateMetaData<double>(meta, "NRRD_thicknesses[2]", spacing);

  std::vector<std::vector<double> > msrFrame(3);
  for( unsigned int saxi = 0; saxi < 3; saxi++ )
    {
    msrFrame[saxi].resize(3);
    for( unsigned int saxj = 0; saxj < 3; saxj++ )
      {
      msrFrame[saxi][saxj] = 0.0;
      }
    }
  msrFrame[0][0] = 1.0; msrFrame[1][1] = 1.0; msrFrame[2][2] = 1.0;
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(meta, "NRRD_measurement frame", msrFrame);

  NrrdValue = "DWMRI";
  itk::EncapsulateMetaData<std::string>(meta, "modality", NrrdValue);

  itk::NumberToString<double> doubleConvert;

  NrrdValue = doubleConvert(m_BValue);
  itk::EncapsulateMetaData<std::string>(meta, "DWMRI_b-value", NrrdValue);
  /* We should apply direction vcl_cosines to gradient directions if requested by
    the user */
  for( int i = 0; i < m_NumberDtiDirections; i++ )
    {
    NrrdValue.clear();

    /* Rotate Diffusion Directions */
    vnl_vector<double> curGradientDirection(3);
    for( int k = 0; k < 3; k++ )
      {
      curGradientDirection[k] = m_DiffusionDirections[i][k];
      }

    if( m_RotateGradients )
      {
      curGradientDirection = m_Input->GetDirection().GetVnlMatrix() * curGradientDirection;
      }
    for( int k = 0; k < 3; k++ )
      {
      NrrdValue += doubleConvert(curGradientDirection[k]);
      if( k < 2 )
        {
        NrrdValue += ' ';
        }
      }
    sprintf(tmpStr, "DWMRI_gradient_%04d", i);
    itk::EncapsulateMetaData<std::string>(meta, tmpStr, NrrdValue);
    }

  m_Output->SetMetaDataDictionary(meta);
}
} // end namespace itk
#endif
