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
//  This class expects the Mesh type to use Scalar as
//  its PixelType (PointData type, to be more specific).
//
//

namespace itk
{
template <class TMesh>
class MeshGeneratorHelper5
{
public:

  typedef TMesh                   MeshType;
  typedef typename TMesh::Pointer MeshPointer;

  static void GenerateMesh( MeshPointer & mesh, unsigned int resolution, double scale, double factor )
  {
    typedef itk::RegularSphereMeshSource<MeshType> SphereMeshSourceType;

    typename SphereMeshSourceType::Pointer sphereMeshSource = SphereMeshSourceType::New();

    typedef typename SphereMeshSourceType::PointType PointType;
    typedef typename PointType::VectorType           VectorType;

    PointType center;
    center.Fill( 0.0 );

    VectorType scaleVector;
    scaleVector.Fill( scale );

    sphereMeshSource->SetCenter( center );
    sphereMeshSource->SetScale( scaleVector );
    sphereMeshSource->SetResolution( resolution );

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

    typedef typename MeshType::PointsContainer     PointContainer;
    typedef typename PointContainer::ConstPointer  PointContainerConstPointer;
    typedef typename PointContainer::ConstIterator PointContainerConstIterator;

    PointDataIterator pixelIterator = pointData->Begin();
    PointDataIterator pixelEnd      = pointData->End();

    PointContainerConstPointer  points = mesh->GetPoints();
    PointContainerConstIterator pointItr = points->Begin();

    while( pixelIterator != pixelEnd )
      {
      const PointType & p = pointItr.Value();

      const double x = p[0];
      const double y = p[1];
      const double z = p[2];

      const double xy = vcl_sqrt( x * x + y * y );

      const double phi = vcl_atan2( xy, z );

      const double Gx = phi * factor;

      pixelIterator.Value() = vcl_exp( -Gx * Gx );

      ++pointItr;
      ++pixelIterator;
      }
  }
};
}

#endif
