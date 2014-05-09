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

#ifndef __itkGtractTransformIO_cxx
#define __itkGtractTransformIO_cxx

#include "itkGtractTransformIO.h"

#include <iostream>

namespace itk
{
GtractTransformIO
::GtractTransformIO()
{
}

void GtractTransformIO::SetFileName( char *fileName)
{
  m_FileName = fileName;
}

void GtractTransformIO::SetFileName( std::string fileName)
{
  m_FileName = fileName;
}

void GtractTransformIO::LoadTransform()
{
  typedef itk::TransformFileReader TransformReaderType;
  TransformReaderType::Pointer transformReader =  TransformReaderType::New();

  transformReader->SetFileName( m_FileName.c_str() );

  std::cout << "Reading ITK transform file: " << m_FileName << " ..." << std::endl;

  try
    {
    transformReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Failed to load Transform: " << std::endl;
    std::cerr << err << std::endl;
    throw;
    }

  std::cout << "Read ITK transform file: " << m_FileName << std::endl;

  std::string readTransformType = ( transformReader->GetTransformList()->back() )->GetTransformTypeAsString();
  // std::cerr << "Transform has type string '" << readTransformType << "'." <<
  // std::endl;

  if( strcmp(readTransformType.c_str(), "VersorRigid3DTransform_double_3_3") == 0 )
    {
    // Load Versor Transform
    m_RigidTransform = RigidTransformType::New();
    m_RigidTransform->SetIdentity();
    m_RigidTransform->SetParameters(
      ( *transformReader->GetTransformList()->begin() )->GetParameters() );
    m_RigidTransform->SetFixedParameters(
      ( *transformReader->GetTransformList()->begin() )->GetFixedParameters() );
    // std::cout << "Parameters " <<
    // (*transformReader->GetTransformList()->begin())->GetParameters() <<
    // std::endl;
    // std::cout << "Fixed Parameters " <<
    // (*transformReader->GetTransformList()->begin())->GetFixedParameters() <<
    // std::endl;
    }
  else if( strcmp(readTransformType.c_str(), "BSplineTransform_double_3_3") == 0 )
    {
    // std::cerr << "Transform with type string '" << readTransformType << "'
    // was recognized." << std::endl;
    // Load B-Spline Transform
    m_BSplineTransform = BSplineTransformType::New();

    m_BSplineTransform->SetFixedParameters(
      ( transformReader->GetTransformList()->back() )->GetFixedParameters() );
    m_BSplineTransform->SetParameters(
      ( transformReader->GetTransformList()->back() )->GetParameters() );
    std::string initTransformType = ( *transformReader->GetTransformList()->begin() )->GetTransformTypeAsString();
    // std::cout << "initTransformType: " << initTransformType << std::endl;
    if( strcmp(initTransformType.c_str(), "VersorRigid3DTransform_double_3_3") == 0 )
      {
      // std::cout << " Set Bulk Transform: " << std::endl;
      RigidTransformType::Pointer bulkTransform = RigidTransformType::New();
      bulkTransform->SetIdentity();
      bulkTransform->SetParameters(
        ( *transformReader->GetTransformList()->begin() )->GetParameters() );
      bulkTransform->SetFixedParameters(
        ( *transformReader->GetTransformList()->begin() )->GetFixedParameters() );
      m_BSplineTransform->SetBulkTransform( bulkTransform );
      }

    // std::cout << " Transform: " << m_BSplineTransform << std::endl;
    }
  else if( strcmp(readTransformType.c_str(), "AffineTransform_double_3_3") == 0 )
    {
    // Import Affine Transform
    m_AffineTransform = AffineTransformType::New();
    m_AffineTransform->SetIdentity();
    m_AffineTransform->SetParameters(
      ( *transformReader->GetTransformList()->begin() )->GetParameters() );
    }
  else if( strcmp(readTransformType.c_str(), "ThinPlateR2LogRSplineKernelTransform_double_3_3") == 0 )
    {
    // Load Thin Plate Spline Transform - B-Spline Inverse
    m_InverseBSplineTransform = ThinPlateSplineTransformType::New();
    std::cout << "Set Fixed Parameters" << std::endl;
    m_InverseBSplineTransform->SetFixedParameters(
      ( *transformReader->GetTransformList()->begin() )->GetFixedParameters() );
    std::cout << "Set Parameters" << std::endl;
    m_InverseBSplineTransform->SetParameters(
      ( *transformReader->GetTransformList()->begin() )->GetParameters() );
    std::cout << "Compute W Matrix" << std::endl;
    m_InverseBSplineTransform->ComputeWMatrix();
    std::cout << "Done" << std::endl;
    }
  else
    {
    std::cerr << "Transform with type string '" << readTransformType << "' was not recognized." << std::endl;
    }
}

void GtractTransformIO::SaveTransform( int type )
{
  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer transformWriter =  TransformWriterType::New();
  transformWriter->SetFileName( m_FileName.c_str() );

  switch( type )
    {
    case 0:
      {
      transformWriter->SetInput( m_RigidTransform );
      }
      break;
    case 1:
      {
      transformWriter->SetInput( m_BSplineTransform->GetBulkTransform() );
      transformWriter->AddTransform( m_BSplineTransform );
      }
      break;
    case 2:
      {
      transformWriter->SetInput( m_AffineTransform );
      }
      break;
    case 3:
      {
      transformWriter->SetInput( m_InverseBSplineTransform );
      }
      break;
    }
  transformWriter->Update();
  std::cout << "Wrote ITK transform to file: " << m_FileName << std::endl;
}
} // end namespace itk
#endif
