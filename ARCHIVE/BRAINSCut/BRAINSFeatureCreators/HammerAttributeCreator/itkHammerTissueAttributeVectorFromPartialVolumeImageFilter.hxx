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

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkHammerTissueAttributeVectorFromPartialVolumeImageFilter.hxx,v $
Language:  C++
Date:      $Date: 2009/01/14 21:46:50 $
Version:   $Revision: 1.6 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for
details.

This program is developed under NIH NCBC collaboration grant
R01 EB006733, "Development and Dissemination of Robust Brain MRI
Measurement Tools".

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHammerTissueAttributeVectorFromPartialVolumeImageFilter_hxx
#define __itkHammerTissueAttributeVectorFromPartialVolumeImageFilter_hxx
#include "itkHammerTissueAttributeVectorFromPartialVolumeImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkProgressReporter.h"
#include "itkIntensityWindowingImageFilter.h"

namespace itk
{
//
// Constructor
//
template <typename TInputImage, typename TOutputImage>
HammerTissueAttributeVectorFromPartialVolumeImageFilter<TInputImage, TOutputImage>
::HammerTissueAttributeVectorFromPartialVolumeImageFilter()
{
  this->m_UseImageSpacing   = true;
  this->m_Scale = 5;
  this->m_Strength = 1;
  this->m_GMValue = 150;
  this->m_WMValue = 250;
  this->m_CSFValue = 10;
  this->m_VNValue = 50;
  this->m_BGValue = 0;

  this->m_OffsetInSphericalNeighborhood.clear();

#if defined( ITK_IMAGE_BEHAVES_AS_ORIENTED_IMAGE )
  this->m_UseImageDirection = true;
#else
  this->m_UseImageDirection = false;
#endif
}

//
// Destructor
//
template <typename TInputImage, typename TOutputImage>
HammerTissueAttributeVectorFromPartialVolumeImageFilter<TInputImage, TOutputImage>
::~HammerTissueAttributeVectorFromPartialVolumeImageFilter()
{
}

template <typename TInputImage, typename TOutputImage>
void
HammerTissueAttributeVectorFromPartialVolumeImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
throw (InvalidRequestedRegionError)
{
  printf("* GenerateInputRequestRegion() \n");
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  InputImagePointer inputPtr =
    const_cast<InputImageType *>( this->GetInput() );
  OutputImagePointer outputPtr = this->GetOutput();

  if( !inputPtr || !outputPtr )
    {
    return;
    }

  // pad input requested region by 1, so the simple edge detection
  // works on the entire image domain
  unsigned long radius = 1;
  typename InputImageType::SpacingType inputSpacing = this->GetInput()->GetSpacing();
  Size<InputImageDimension> sphereRadius;
  for( unsigned int i = 0; i < InputImageDimension; ++i )
    {
    sphereRadius[i] = static_cast<unsigned long>( this->m_Scale / inputSpacing[i] );
    if( sphereRadius[i] > radius )
      {
      radius = sphereRadius[i];
      }
    }

  // compute spherical neighborhood for geometrical attribute
  // computation
  typename InputImageType::Pointer dummyImage = InputImageType::New();
  typename InputImageType::RegionType dummyRegion;
  typename InputImageType::PointType dummyOrigin;
  typename InputImageType::IndexType dummyStart;
  typename InputImageType::SizeType dummySize;

  dummyImage->SetSpacing(inputSpacing);
  for( unsigned int k = 0; k < InputImageDimension; k++ )
    {
    dummySize[k] = sphereRadius[k] + sphereRadius[k] + 1;
    dummyStart[k] = -sphereRadius[k];
    dummyOrigin[k] = 0;
    }
  dummyRegion.SetIndex(dummyStart);
  dummyRegion.SetSize(dummySize);
  dummyImage->SetRegions(dummyRegion);
  dummyImage->SetOrigin(dummyOrigin);

  float                                             radiusSqr = this->m_Scale * this->m_Scale;
  itk::ImageRegionIteratorWithIndex<InputImageType> it(dummyImage, dummyRegion);
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    typename InputImageType::IndexType dummyIdx = it.GetIndex();
    typename InputImageType::PointType point;
    dummyImage->TransformIndexToPhysicalPoint(dummyIdx, point);
    float d = 0;
    for( unsigned int k = 0; k < InputImageDimension; k++ )
      {
      d += point[k] * point[k];
      }
    if( d > radiusSqr )
      {
      continue;
      }
    else
      {
      NeighborOffsetType offset;
      for( unsigned int k = 0; k < InputImageDimension; k++ )
        {
        offset[k] = dummyIdx[k];
        }
      this->m_OffsetInSphericalNeighborhood.push_back(offset);
      }
    }

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius(radius);

  // crop the input requested region at the input's largest possible region
  if( inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() ) )
    {
    inputPtr->SetRequestedRegion(inputRequestedRegion);
    return;
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion(inputRequestedRegion);

    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);

    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
    }
}

template <typename TInputImage, typename TOutputImage>
void
HammerTissueAttributeVectorFromPartialVolumeImageFilter<TInputImage, TOutputImage>
::CreateN1Neighbor()
{
  printf("* CreateN1Neighbor() \n");
  m_N1Neighborhood.resize(6);
  for( int k = 0; k < 6; k++ )
    {
    for( size_t s = 0; s < InputImageDimension; s++ )
      {
      m_N1Neighborhood[k][s] = 0;
      }
    }
  m_N1Neighborhood[0][0] = -1;
  m_N1Neighborhood[1][0] = 1;
  m_N1Neighborhood[2][1] = -1;
  m_N1Neighborhood[3][1] = 1;
  m_N1Neighborhood[4][2] = -1;
  m_N1Neighborhood[5][2] = 1;
}

template <typename TInputImage, typename TOutputImage>
void
HammerTissueAttributeVectorFromPartialVolumeImageFilter<TInputImage, TOutputImage>
::CreateFeatureNeighbor(int Radius)
{
  printf("* CreateFeatureNeighbor\n");
  // InputSpacingType spacing = this->GetInput()->GetSpacing();
  float              rad_sqrd = Radius * Radius;
  int                i, j, k;
  NeighborOffsetType offset;

  m_FeatureNeighborhood.clear();
  for( i = -Radius; i <= Radius; i++ )
    {
    for( j = -Radius; j <= Radius; j++ )
      {
      for( k = -Radius; k <= Radius; k++ )
        {
        if( ( i * i + j * j + k * k ) <= rad_sqrd )
          {
          offset[0] = static_cast<long int>( i );
          offset[1] = static_cast<long int>( j );
          offset[2] = static_cast<long int>( float(k) / 1.5 );
          m_FeatureNeighborhood.push_back(offset);
          }
        }
      }
    }
}

template <typename TInputImage, typename TOutputImage>
void
HammerTissueAttributeVectorFromPartialVolumeImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  printf("* GenerateData() \n ");
  this->AllocateOutputs();

  typename TOutputImage::PixelType attributeVector;

  // Get the input and output
  OutputImageType * outputImage = this->GetOutput();

  const InputImageType *inputGMVolume = this->GetInput(0);
  const InputImageType *inputWMVolume = this->GetInput(1);
  const InputImageType *inputCSFVolume = this->GetInput(2);

  // Use inputVolume as a reference volume for image information
  //

  InputRegionType dummyRegion = inputWMVolume->GetLargestPossibleRegion();
  // Create the neighbor
  CreateN1Neighbor();
  CreateFeatureNeighbor( static_cast<int>( m_Scale ) );

  // Initialize
  outputImage->SetLargestPossibleRegion(
    inputWMVolume->GetLargestPossibleRegion() );
  ImageRegionIteratorWithIndex<OutputImageType> it( outputImage, outputImage->GetLargestPossibleRegion() );

  ImageRegionConstIteratorWithIndex<InputImageType> source( inputWMVolume, inputWMVolume->GetLargestPossibleRegion() );
  attributeVector.Fill( 0 );
  for( it.GoToBegin(), source.GoToBegin(); !it.IsAtEnd(); ++it, ++source )
    {
    //
    // copy the WM Posterior to the attirubuteVector[1]
    // For WHAT though??????
    //
    attributeVector[1] = source.Get();
    it.Set(attributeVector);
    }

  // Compute the Edge Information
  float flag_GM, flag_CSF, flag_WM;
  //
  // TODO make strength as a input parameter.
  //
  int strength = 1;
  typename InputImageType::PixelType centerPixel;
  typename InputImageType::PixelType PWm, PGm, PCsf;
  std::string    centerTissueType = ""; // WM, GM, or CSF
  InputIndexType centerIdx, neighborIdx;
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    //
    // [WM edge]
    // For voxle WM > 0.3, in 6 neghibors
    //     if #( GM > 0.3) > #( CSF > 0.3) ==> WM || GM edge
    //     else if #( CSF > 0.3) > 0        ==> WM || CSF edge
    //
    attributeVector = it.Get();
    centerIdx = it.GetIndex();

    // determin which tissue type it belongs to.
    //

    PWm = inputWMVolume->GetPixel( centerIdx );
    PGm = inputGMVolume->GetPixel( centerIdx );
    PCsf = inputCSFVolume->GetPixel( centerIdx );

    if( PWm > PGm )
      {
      if( PWm > PCsf )
        {
        centerTissueType = "WM"; centerPixel = PWm;
        }
      else
        {
        centerTissueType = "CSF"; centerPixel = PCsf;
        }
      }
    else
      {
      if( PGm > PCsf )
        {
        centerTissueType = "GM"; centerPixel = PGm;
        }
      else
        {
        centerTissueType = "CSF"; centerPixel = PCsf;
        }
      }

    if( centerPixel < 0.1 )
      {
      centerTissueType = "NONE";
      centerPixel = 0;
      }

    // Determin which Bondary the edge is
    //
    if( centerIdx[0] == 0 || centerIdx[1] == 0 || centerIdx[2] == 0 )
      {
      continue;
      }
    if( centerIdx[0] == static_cast<signed int>(dummyRegion.GetSize()[0] - 1) ||
        centerIdx[1] == static_cast<signed int>(dummyRegion.GetSize()[1] - 1) ||
        centerIdx[2] == static_cast<signed int>(dummyRegion.GetSize()[2] - 1) )
      {
      continue;
      }
    flag_GM = 0.0F;
    flag_CSF = 0.0F;
    flag_WM = 0.0F;
    for( unsigned int k = 0; k < m_N1Neighborhood.size(); k++ )
      {
      for( size_t s = 0; s < InputImageDimension; s++ )
        {
        neighborIdx[s] = centerIdx[s] + m_N1Neighborhood[k][s];
        }

      // Sum
      flag_GM += inputGMVolume->GetPixel( neighborIdx);
      flag_CSF += inputCSFVolume->GetPixel( neighborIdx);
      flag_WM += inputWMVolume->GetPixel( neighborIdx);
      }
    if( centerTissueType == "WM" )
      {
      if( flag_GM > flag_CSF && flag_GM > strength )
        {
        attributeVector[0] = m_WMGMEDGE;
        it.Set(attributeVector);
        }
      if( flag_GM <= flag_CSF && flag_CSF > strength )
        {
        attributeVector[0] = m_WMCSFEDGE;
        it.Set(attributeVector);
        }
      }
    if( centerTissueType == "GM" )
      {
      if( flag_WM > flag_CSF && flag_WM > strength )
        {
        attributeVector[0] = m_WMGMEDGE;
        it.Set(attributeVector);
        }
      if( flag_WM <= flag_CSF && flag_CSF > strength )
        {
        attributeVector[0] = m_WMCSFEDGE;
        it.Set(attributeVector);
        }
      }
    if( centerTissueType == "CSF" )
      {
      if( flag_GM > flag_WM && flag_GM > strength )
        {
        attributeVector[0] = m_WMCSFEDGE;
        it.Set(attributeVector);
        }
      if( flag_GM <= flag_WM && flag_GM > strength )
        {
        attributeVector[0] = m_GMCSFEDGE;
        it.Set(attributeVector);
        }
      }
    } // end of Edge computation

  //
  // [Compute the GMIs]
  //
  float NonWM_value, CSF_value, GM_value;
  float pixelNumInBubble = m_FeatureNeighborhood.size();
  printf("* Here we compute GMIs\n");
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    attributeVector = it.Get();
    centerIdx = it.GetIndex();

    // float sum = inputWMVolume->GetPixel(centerIdx) +
    //            inputGMVolume->GetPixel(centerIdx) +
    //            inputCSFVolume->GetPixel(centerIdx) ;
    // if( sum < 0.1F )
    if( attributeVector.GetEdge() == 0 )
      {
      continue;
      }
    NonWM_value = 0.0F;
    CSF_value = 0.0F;
    GM_value = 0.0F;
    for( unsigned int t = 0; t < m_FeatureNeighborhood.size(); t++ )
      {
      for( size_t s = 0; s < InputImageDimension; s++ )
        {
        neighborIdx[s] = centerIdx[s] + (int)m_FeatureNeighborhood[t][s];
        }
      if( neighborIdx[0] < 0 || neighborIdx[1] < 0 || neighborIdx[2] < 0 )
        {
        continue;
        }

      if( neighborIdx[0] >= static_cast<signed int>(dummyRegion.GetSize()[0]) ||
          neighborIdx[1] >= static_cast<signed int>(dummyRegion.GetSize()[1]) ||
          neighborIdx[2] >= static_cast<signed int>(dummyRegion.GetSize()[2]) )
        {
        continue;
        }

      // Sum all the probaiblity in the ball
      //
      NonWM_value += (1.0F - inputWMVolume->GetPixel(neighborIdx) );
      CSF_value   += (float)(inputCSFVolume->GetPixel(neighborIdx) );
      GM_value    += (float)(inputGMVolume->GetPixel(neighborIdx) );

      float degree = (NonWM_value / pixelNumInBubble);
      attributeVector[2] = degree * 100.0F;

      float CSF_degree = CSF_value / pixelNumInBubble;
      attributeVector[3] = CSF_degree * 100.0F;

      float GM_degree = GM_value / pixelNumInBubble;
      attributeVector[4] = GM_degree * 100.0F;

      it.Set(attributeVector);
      }
    }
}

/**
  * Standard "PrintSelf" method
  */
template <typename TInputImage, typename TOutputImage>
void
HammerTissueAttributeVectorFromPartialVolumeImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "UseImageSpacing: "
     << ( this->m_UseImageSpacing ? "On" : "Off" ) << std::endl;
  os << indent << "UseImageDirection = "
     << ( this->m_UseImageDirection ? "On" : "Off" ) << std::endl;
}
}   // end namespace itk

#endif
