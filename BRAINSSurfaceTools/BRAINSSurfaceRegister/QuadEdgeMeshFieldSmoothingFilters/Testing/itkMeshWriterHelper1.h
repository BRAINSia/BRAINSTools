/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeanSquaresMeshToMeshMetricTest1.cxx,v $
  Language:  C++
  Date:      $Date: 2007-09-06 17:44:24 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkMeshWriterHelper1_h
#define __itkMeshWriterHelper1_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

namespace itk
{
template <class TMesh>
class MeshWriterHelper1
{
public:
  static void WriteMeshToFile( const TMesh * mesh, const char * filename )
  {
    typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<TMesh> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( mesh );
    writer->SetFileName( filename );

    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Error during writer Update() " << std::endl;
      std::cerr << excp << std::endl;
      return;
      }
  }
};
}

#endif
