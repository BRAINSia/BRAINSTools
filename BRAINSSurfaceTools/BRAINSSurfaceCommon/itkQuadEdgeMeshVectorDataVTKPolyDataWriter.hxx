/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkQuadEdgeMeshVectorDataVTKPolyDataWriter.txx,v $
  Language:  C++
  Date:      $Date: 2009-01-02 18:43:05 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkQuadEdgeMeshVectorDataVTKPolyDataWriter_hxx
#define __itkQuadEdgeMeshVectorDataVTKPolyDataWriter_hxx

#include "itkQuadEdgeMeshVectorDataVTKPolyDataWriter.h"
#include "itkMeasurementVectorTraits.h"
#include "itkNumberToString.h"

namespace itk
{
//
// Constructor
//
template <class TMesh>
QuadEdgeMeshVectorDataVTKPolyDataWriter<TMesh>
::QuadEdgeMeshVectorDataVTKPolyDataWriter()
{
  m_CellDataName = "";
  m_PointDataName = "";
}

//
// Destructor
//
template <class TMesh>
QuadEdgeMeshVectorDataVTKPolyDataWriter<TMesh>
::~QuadEdgeMeshVectorDataVTKPolyDataWriter()
{
}

template <class TMesh>
void
QuadEdgeMeshVectorDataVTKPolyDataWriter<TMesh>
::GenerateData()
{
  this->Superclass::GenerateData();
  this->WriteCellData();
  this->WritePointData();
}

template <class TMesh>
void
QuadEdgeMeshVectorDataVTKPolyDataWriter<TMesh>
::WriteCellData()
{
  CellDataContainerConstPointer celldata = this->m_Input->GetCellData();
  itk::NumberToString<double>                doubleToString;

  if( celldata )
    {
    if( celldata->Size() != 0 )
      {
      std::ofstream outputFile( this->m_FileName.c_str(), std::ios_base::app );

      outputFile << "CELL_DATA " << this->m_Input->GetNumberOfFaces() << std::endl;
      outputFile << "SCALARS ";

      if( m_CellDataName != "" )
        {
        outputFile << m_CellDataName << " " << m_CellDataName << std::endl;
        }
      else
        {
        outputFile << "double double" << std::endl;
        }

      outputFile << "LOOKUP_TABLE default" << std::endl;

      unsigned long k(0);

      CellsContainerConstPointer  cells = this->m_Input->GetCells();
      CellsContainerConstIterator it = cells->Begin();

      CellDataContainerConstIterator c_it = celldata->Begin();

      while( c_it != celldata->End() )
        {
        CellType* cellPointer = it.Value();
        if( cellPointer->GetType() != 1 )
          {
          for( unsigned _i = 0; _i < c_it.Value().Size(); ++_i )
            {
            outputFile << doubleToString(c_it.Value()[_i]);
            if( _i < 2 )
              {
              outputFile << " ";
              }
            }
          if( k++ % 3 == 0 )
            {
            outputFile << std::endl;
            }
          }
        ++c_it;
        ++it;
        }

      outputFile << std::endl;
      outputFile.close();
      }
    }
}

template <class TMesh>
void
QuadEdgeMeshVectorDataVTKPolyDataWriter<TMesh>
::WritePointData()
{
  PointDataContainerConstPointer pointdata = this->m_Input->GetPointData();
  itk::NumberToString<double>                 doubleToString;

  if( pointdata )
    {
    std::ofstream outputFile( this->m_FileName.c_str(), std::ios_base::app );

    outputFile << "POINT_DATA " << this->m_Input->GetNumberOfPoints() << std::endl;
    outputFile << "VECTORS ";

    if( m_PointDataName != "" )
      {
      outputFile << m_PointDataName << " double" << std::endl;
      }
    else
      {
      outputFile << "vectors double" << std::endl;
      }

    unsigned long k = 0;

    PointDataContainerIterator c_it = pointdata->Begin();

    const unsigned int vectorSize = c_it.Value().Size(); // Statistics::MeasurementVectorTraits::GetLength( c_it.Value()
                                                         // );

    while(  c_it != pointdata->End() )
      {
      for( unsigned int j = 0; j < vectorSize; j++ )
        {
        outputFile << doubleToString(static_cast<double>( c_it.Value()[j] ) ) << " ";
        ++k;
        if( k % 3 == 0 )
          {
          outputFile << std::endl;
          }
        }

      ++c_it;
      }

    outputFile << std::endl;
    outputFile.close();
    }
}
}

#endif
