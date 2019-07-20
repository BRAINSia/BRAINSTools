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

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkInvertBSplineFilter.h"
#include <itkIOCommon.h>
#include <itkCastImageFilter.h>
#include <itkConstantPadImageFilter.h>
#include <itkExtractImageFilter.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"

#include <iostream>

namespace itk
{
InvertBSplineFilter::InvertBSplineFilter()
{
  m_XgridSize = 2;
  m_YgridSize = 2;
  m_ZgridSize = 2;
}

void
InvertBSplineFilter::Update()
{
  std::cout << "InvertBSplineFilter()...." << std::endl;

  ImageType::SizeType   imageSize = m_ExampleImage->GetLargestPossibleRegion().GetSize();
  PointSetType::Pointer sourceLandMarks = PointSetType::New();
  PointSetType::Pointer targetLandMarks = PointSetType::New();

  PointSetType::PointsContainer::Pointer sourceLandMarkContainer = sourceLandMarks->GetPoints();
  PointSetType::PointsContainer::Pointer targetLandMarkContainer = targetLandMarks->GetPoints();

  float xinr = (float)( imageSize[0] ) / (float)m_XgridSize;
  float yinr = (float)( imageSize[1] ) / (float)m_YgridSize;
  float zinr = (float)( imageSize[2] ) / (float)m_ZgridSize;

  PointIdType id = itk::NumericTraits< PointIdType >::ZeroValue();
  std::cout << "Xsize " << imageSize[0] << " Ysize " << imageSize[1] << " Zsize " << imageSize[2] << std::endl;
  std::cout << "Xinc " << xinr << " Yinc " << yinr << " Zinc " << zinr << std::endl;

  for ( float z = 0; z <= imageSize[2] + 1; z += zinr )
  {
    for ( float y = 0; y <= imageSize[1] + 1; y += yinr )
    {
      for ( float x = 0; x <= imageSize[0] + 1; x += xinr )
      {
        PointType                         p1;
        PointType                         p2;
        itk::ContinuousIndex< double, 3 > imageIndex;
        imageIndex[0] = x;
        imageIndex[1] = y;
        imageIndex[2] = z;
        m_ExampleImage->TransformContinuousIndexToPhysicalPoint( imageIndex, p1 );
        p2 = m_Input->TransformPoint( p1 );

        // Add Points and the Landmark Point Lists
        sourceLandMarkContainer->InsertElement( id, p2 );
        targetLandMarkContainer->InsertElement( id, p1 );
        id++;
        // std::cout << "Set Z points " << z << " " << p1 << " " << p2 <<
        // std::endl;
      }
    }
  }
  std::cout << "Set Landmarks and Invert " << std::endl;

  // Create TPS Approximation of B-Spline Transform
  m_Output = TransformType::New();
  m_Output->SetSourceLandmarks( sourceLandMarks );
  m_Output->SetTargetLandmarks( targetLandMarks );
  m_Output->ComputeWMatrix();

  std::cout << "Computed TPS Inverse transform" << std::endl;

  std::cout << m_Output << std::endl;
}
} // end namespace itk
