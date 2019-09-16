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
Module:    $RCSfile: itkHammerTissueAttributeVectorImageFilter.hxx,v $
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
#ifndef __itkHammerTissueAttributeVectorImageFilter_hxx
#define __itkHammerTissueAttributeVectorImageFilter_hxx
#include "itkHammerTissueAttributeVectorImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkProgressReporter.h"

namespace itk
{
//
// Constructor
//
template <typename TInputImage, typename TOutputImage>
HammerTissueAttributeVectorImageFilter<TInputImage, TOutputImage>::HammerTissueAttributeVectorImageFilter()
{
  this->m_UseImageSpacing = true;
  this->m_Scale = 5;
  this->m_Strength = 1;
  this->m_GMValue = 150;
  this->m_WMValue = 250;
  this->m_CSFValue = 10;
  this->m_VNValue = 50;
  this->m_BGValue = 0;

  this->m_OffsetInSphericalNeighborhood.clear();

#if defined(ITK_IMAGE_BEHAVES_AS_ORIENTED_IMAGE)
  this->m_UseImageDirection = true;
#else
  this->m_UseImageDirection = false;
#endif
}

//
// Destructor
//
template <typename TInputImage, typename TOutputImage>
HammerTissueAttributeVectorImageFilter<TInputImage, TOutputImage>::~HammerTissueAttributeVectorImageFilter()
{}

template <typename TInputImage, typename TOutputImage>
void
HammerTissueAttributeVectorImageFilter<TInputImage, TOutputImage>::GenerateInputRequestedRegion() throw(
  InvalidRequestedRegionError)
{
  printf("* GenerateInputRequestRegion() \n");
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  InputImagePointer  inputPtr = const_cast<InputImageType *>(this->GetInput());
  OutputImagePointer outputPtr = this->GetOutput();

  if (!inputPtr || !outputPtr)
  {
    return;
  }

  // pad input requested region by 1, so the simple edge detection
  // works on the entire image domain
  unsigned long                        radius = 1;
  typename InputImageType::SpacingType inputSpacing = this->GetInput()->GetSpacing();
  Size<InputImageDimension>            sphereRadius;
  for (int i = 0; i < InputImageDimension; ++i)
  {
    sphereRadius[i] = static_cast<unsigned long>(this->m_Scale / inputSpacing[i]);
    if (sphereRadius[i] > radius)
    {
      radius = sphereRadius[i];
    }
  }

  // compute spherical neighborhood for geometrical attribute
  // computation
  typename InputImageType::Pointer    dummyImage = InputImageType::New();
  typename InputImageType::RegionType dummyRegion;
  typename InputImageType::PointType  dummyOrigin;
  typename InputImageType::IndexType  dummyStart;
  typename InputImageType::SizeType   dummySize;

  dummyImage->SetSpacing(inputSpacing);
  for (int k = 0; k < InputImageDimension; k++)
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
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    typename InputImageType::IndexType dummyIdx = it.GetIndex();
    typename InputImageType::PointType point;
    dummyImage->TransformIndexToPhysicalPoint(dummyIdx, point);
    float d = 0;
    for (int k = 0; k < InputImageDimension; k++)
    {
      d += point[k] * point[k];
    }
    if (d > radiusSqr)
    {
      continue;
    }
    else
    {
      NeighborOffsetType offset;
      for (int k = 0; k < InputImageDimension; k++)
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
  if (inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()))
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
HammerTissueAttributeVectorImageFilter<TInputImage, TOutputImage>::CreateN1Neighbor()
{
  printf("* CreateN1Neighbor() \n");
  m_N1Neighborhood.resize(6);
  for (int k = 0; k < 6; k++)
  {
    for (int s = 0; s < InputImageDimension; s++)
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
HammerTissueAttributeVectorImageFilter<TInputImage, TOutputImage>::CreateFeatureNeighbor(int Radius)
{
  printf("* CreateFeatureNeighbor\n");
  // InputSpacingType spacing = this->GetInput()->GetSpacing();
  float              rad_sqrd = Radius * Radius;
  int                i, j, k;
  NeighborOffsetType offset;

  m_FeatureNeighborhood.clear();
  for (i = -Radius; i <= Radius; i++)
  {
    for (j = -Radius; j <= Radius; j++)
    {
      for (k = -Radius; k <= Radius; k++)
      {
        if ((i * i + j * j + k * k) <= rad_sqrd)
        {
          offset[0] = static_cast<long int>(i);
          offset[1] = static_cast<long int>(j);
          offset[2] = static_cast<long int>(float(k) / 1.5);
          m_FeatureNeighborhood.push_back(offset);
        }
      }
    }
  }
}

template <typename TInputImage, typename TOutputImage>
void
HammerTissueAttributeVectorImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  printf("* GenerateData() \n ");
  this->AllocateOutputs();

  typename TOutputImage::PixelType attributeVector;

  // Get the input and output
  OutputImageType *                    outputImage = this->GetOutput();
  const InputImageType *               inputImage = this->GetInput();
  typename InputImageType::SpacingType inputSpacing = inputImage->GetSpacing();
  InputRegionType                      dummyRegion = inputImage->GetLargestPossibleRegion();
  // Create the neighbor
  CreateN1Neighbor();
  CreateFeatureNeighbor(static_cast<int>(m_Scale));

  // Initialize
  ImageRegionIteratorWithIndex<OutputImageType> it(outputImage, outputImage->GetLargestPossibleRegion());

  ImageRegionConstIteratorWithIndex<InputImageType> source(inputImage, inputImage->GetLargestPossibleRegion());
  attributeVector.Fill(0);
  int voxel_num = 0;
  for (it.GoToBegin(), source.GoToBegin(); !it.IsAtEnd(); ++it, ++source)
  {
    attributeVector[1] = source.Get();
    it.Set(attributeVector);
    if (attributeVector[1] == m_WMValue)
    {
      voxel_num++;
    }
  }
  printf(" voxel_num=%d\n", voxel_num);

  // Compute the Edge Information
  int                                flag_GM, flag_VN, flag_CSF;
  int                                strength = 1;
  int                                edge_num = 0;
  typename InputImageType::PixelType centerPixel, currentPixel;
  InputIndexType                     idx, cur_idx;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    attributeVector = it.Get();
    centerPixel = attributeVector[1];
    idx = it.GetIndex();
    if (idx[0] == 0 || idx[1] == 0 || idx[2] == 0)
    {
      continue;
    }
    if (idx[0] == static_cast<signed int>(dummyRegion.GetSize()[0] - 1) ||
        idx[1] == static_cast<signed int>(dummyRegion.GetSize()[1] - 1) ||
        idx[2] == static_cast<signed int>(dummyRegion.GetSize()[2] - 1))
    {
      continue;
    }
    if (centerPixel == m_WMValue)
    {
      flag_GM = 0;
      flag_VN = 0;
      for (unsigned int k = 0; k < m_N1Neighborhood.size(); k++)
      {
        for (int s = 0; s < InputImageDimension; s++)
        {
          cur_idx[s] = idx[s] + m_N1Neighborhood[k][s];
        }

        //        if(!dummyRegion.IsInsideInWorldSpace(cur_idx))
        //          continue;

        currentPixel = inputImage->GetPixel(cur_idx);
        if (currentPixel == m_GMValue)
        {
          flag_GM++;
        }
        if (currentPixel == m_VNValue)
        {
          flag_VN++;
        }
      }

      /* determine, whether edge? If edge, which type */
      if (flag_GM > flag_VN)
      {
        attributeVector[0] = m_WMGMEDGE;
        it.Set(attributeVector);
      }
      if (flag_GM <= flag_VN && flag_VN != 0)
      {
        attributeVector[0] = m_WMVNEDGE;
        it.Set(attributeVector);
      }
    }
  }
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    attributeVector = it.Get();
    centerPixel = attributeVector[1];
    idx = it.GetIndex();
    if (idx[0] == 0 || idx[1] == 0 || idx[2] == 0)
    {
      continue;
    }
    if (idx[0] == static_cast<signed int>(dummyRegion.GetSize()[0] - 1) ||
        idx[1] == static_cast<signed int>(dummyRegion.GetSize()[1] - 1) ||
        idx[2] == static_cast<signed int>(dummyRegion.GetSize()[2] - 1))
    {
      continue;
    }
    if (centerPixel == m_GMValue)
    {
      flag_CSF = 0;
      for (unsigned int k = 0; k < m_N1Neighborhood.size(); k++)
      {
        for (int s = 0; s < InputImageDimension; s++)
        {
          cur_idx[s] = idx[s] + m_N1Neighborhood[k][s];
        }
        //        if(!dummyRegion.IsInsideInWorldSpace(cur_idx))
        //          continue;
        currentPixel = inputImage->GetPixel(cur_idx);
        if (currentPixel <= m_CSFValue)
        {
          flag_CSF++;
        }
      }
      if (flag_CSF >= strength)
      {
        attributeVector[0] = m_GMCSFBGEDGE;
        it.Set(attributeVector);
      }
    }
  } // end of Edge computation
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    attributeVector = it.Get();
    centerPixel = attributeVector[1];
    idx = it.GetIndex();
    if (idx[0] == 0 || idx[1] == 0 || idx[2] == 0)
    {
      continue;
    }
    if (idx[0] == static_cast<signed int>(dummyRegion.GetSize()[0] - 1) ||
        idx[1] == static_cast<signed int>(dummyRegion.GetSize()[1] - 1) ||
        idx[2] == static_cast<signed int>(dummyRegion.GetSize()[2] - 1))
    {
      continue;
    }
    if (centerPixel == m_GMValue)
    {
      flag_VN = 0;
      for (unsigned int k = 0; k < m_N1Neighborhood.size(); k++)
      {
        for (int s = 0; s < InputImageDimension; s++)
        {
          cur_idx[s] = idx[s] + m_N1Neighborhood[k][s];
        }
        //        if(!dummyRegion.IsInsideInWorldSpace(cur_idx))
        //          continue;
        currentPixel = inputImage->GetPixel(cur_idx);
        if (currentPixel == m_VNValue)
        {
          flag_VN++;
        }
      }
      if (flag_VN >= strength)
      {
        attributeVector[0] = m_GMVNEDGE;
        it.Set(attributeVector);
      }
    }
  } // end of Edge computation

  edge_num = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    attributeVector = it.Get();
    if (attributeVector[0] == m_GMCSFBGEDGE)
    {
      edge_num++;
    }
  }
  printf(" edge_num=%d\n", edge_num);
  // Compute the GMIs
  int   NumDown, VN_value, CSFBG_value;
  float pixelNumInBubble = m_FeatureNeighborhood.size();
  printf("* Here we compute GMIs\n");
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    attributeVector = it.Get();
    idx = it.GetIndex();
    if (attributeVector.GetEdge() == 0)
    {
      continue;
    }
    NumDown = 0;
    VN_value = 0;
    CSFBG_value = 0;
    for (unsigned int t = 0; t < m_FeatureNeighborhood.size(); t++)
    {
      for (int s = 0; s < InputImageDimension; s++)
      {
        cur_idx[s] = idx[s] + (int)m_FeatureNeighborhood[t][s];
      }
      //       if(!dummyRegion.IsInsideInWorldSpace(cur_idx))
      //         continue;
      if (cur_idx[0] < 0 || cur_idx[1] < 0 || cur_idx[2] < 0)
      {
        continue;
      }
      if (cur_idx[0] >= static_cast<signed int>(dummyRegion.GetSize()[0]) ||
          cur_idx[1] >= static_cast<signed int>(dummyRegion.GetSize()[1]) ||
          cur_idx[2] >= static_cast<signed int>(dummyRegion.GetSize()[2]))
      {
        continue;
      }
      currentPixel = inputImage->GetPixel(cur_idx);
      if (currentPixel != m_WMValue)
      {
        NumDown++;
      }
      if (currentPixel == m_VNValue)
      {
        VN_value++;
      }
      if (currentPixel <= m_CSFValue)
      {
        CSFBG_value++;
      }
    }
    float degree = ((float)NumDown / pixelNumInBubble);
    float value = degree * 255 * 1.2; /*degree*degree*255*/
    if (value > 255)
    {
      value = 255;
    }
    attributeVector[2] = (unsigned char)value;
    /*printf("(%d,%d,%d)=%d\n", k, i, j, NumDown) ;*/

    /* VN volume */
    VN_value = static_cast<int>(VN_value * 255 * 1.2 / pixelNumInBubble);
    if (VN_value > 255)
    {
      VN_value = 255;
    }
    attributeVector[3] = static_cast<unsigned char>(VN_value);

    /* CSF/BG volume */
    CSFBG_value = static_cast<int>(CSFBG_value * 255 * 1.2 / pixelNumInBubble);
    if (CSFBG_value > 255)
    {
      CSFBG_value = 255;
    }
    attributeVector[4] = static_cast<unsigned char>(CSFBG_value);
    it.Set(attributeVector);
  }
}

/**
 * Standard "PrintSelf" method
 */
template <typename TInputImage, typename TOutputImage>
void
HammerTissueAttributeVectorImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "UseImageSpacing: " << (this->m_UseImageSpacing ? "On" : "Off") << std::endl;
  os << indent << "UseImageDirection = " << (this->m_UseImageDirection ? "On" : "Off") << std::endl;
}
} // end namespace itk

#endif
