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

#ifndef __itkMeshGeneratorHelper2_h
#define __itkMeshGeneratorHelper2_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkQuadEdgeMesh.h"
#include "itkRegularSphereMeshSource.h"

//
//  This class expects the Mesh type to use Vectors as
//  its PixelType (PointData type, to be more specific).
//
//

namespace itk
{
template <class TMesh>
class MeshGeneratorHelper2
{
public:

  typedef TMesh                   MeshType;
  typedef typename TMesh::Pointer MeshPointer;

  static void GenerateMesh( MeshPointer & mesh )
  {
    typedef itk::RegularSphereMeshSource<MeshType> SphereMeshSourceType;

    typename SphereMeshSourceType::Pointer sphereMeshSource = SphereMeshSourceType::New();

    typedef typename SphereMeshSourceType::PointType PointType;
    typedef typename PointType::VectorType           VectorType;

    PointType center;
    center.Fill( 0.0 );

    VectorType scale;
    scale.Fill( 1.0 );

    sphereMeshSource->SetCenter( center );
    sphereMeshSource->SetScale( scale );
    sphereMeshSource->SetResolution( 1 );

    try
      {
      sphereMeshSource->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Error during source Update() " << std::endl;
      std::cerr << excp << std::endl;
      return;
      }

    mesh = sphereMeshSource->GetOutput();

    typedef typename MeshType::PointDataContainer PointDataContainer;

    typename PointDataContainer::Pointer pointData = PointDataContainer::New();

    typedef typename PointDataContainer::Iterator PointDataIterator;

    pointData->Reserve( mesh->GetNumberOfPoints() );

    mesh->SetPointData( pointData );

    PointDataIterator pixelIterator = pointData->Begin();
    PointDataIterator pixelEnd      = pointData->End();

    typedef typename MeshType::PointsContainer::Iterator PointIterator;
    PointIterator pointIterator = mesh->GetPoints()->Begin();
    PointIterator pointEnd      = mesh->GetPoints()->End();

    VectorType vector;
    vector.Fill( 0.0 );

    while( pixelIterator != pixelEnd  && pointIterator != pointEnd )
      {
      const PointType & point = pointIterator.Value();
      const double      x = point[0];
      const double      y = point[1];
      const double      z = point[2];
      const double      w = vcl_sqrt( x * x + y * y );
      const double      phi   = atan2( y, x );
      const double      theta = atan2( z, w );
      vector[0] = vcl_sin( 2 * phi );
      vector[1] = vcl_sin( theta );
      pixelIterator.Value() = vector;
      ++pixelIterator;
      ++pointIterator;
      }
  }
};
}

#endif
