/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFreeSurferBinarySurfaceReader.txx,v $
  Language:  C++
  Date:      $Date: 2008-06-15 02:42:12 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFreeSurferBinarySurfaceReader_txx
#define __itkFreeSurferBinarySurfaceReader_txx

#include "itkFreeSurferBinarySurfaceReader.h"
#include "itkByteSwapper.h"

#include "itksys/SystemTools.hxx"

namespace itk
{
//
// Constructor
//
template <class TOutputMesh>
FreeSurferBinarySurfaceReader<TOutputMesh>
::FreeSurferBinarySurfaceReader()
{
  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, output.GetPointer() );

  this->m_InputGeometryFile = NULL;
  this->m_InputDataFile = NULL;
}

//
// Destructor
//
template <class TOutputMesh>
FreeSurferBinarySurfaceReader<TOutputMesh>
::~FreeSurferBinarySurfaceReader()
{
  this->ReleaseResources();
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReleaseResources()
{
  if( this->m_InputGeometryFile )
    {
    delete this->m_InputGeometryFile;
    }

  if( this->m_InputDataFile )
    {
    delete this->m_InputDataFile;
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::GenerateData()
{
  this->OpenGeometryFile();
  this->ReadGeometryHeader();
  this->PrepareOutputMesh();
  this->ReadSurface();
  this->CloseGeometryFile();

  if( this->DataFileIsAvailable() )
    {
    this->OpenDataFile();
    this->ReadDataHeader();
    this->ReadPointData();
    this->CloseDataFile();
    }
}

template <class TOutputMesh>
bool
FreeSurferBinarySurfaceReader<TOutputMesh>
::DataFileIsAvailable() const
{
  if( this->m_DataFileName.empty() )
    {
    return false;
    }

  if( !itksys::SystemTools::FileExists( m_DataFileName.c_str() ) )
    {
    itkWarningMacro("File " << this->m_DataFileName << " does not exist");
    return false;
    }

  return true;
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::OpenGeometryFile()
{
  if( this->m_FileName.empty() )
    {
    this->ReleaseResources();
    itkExceptionMacro("No input FileName");
    }

  if( !itksys::SystemTools::FileExists( m_FileName.c_str() ) )
    {
    this->ReleaseResources();
    itkExceptionMacro("File " << this->m_FileName << " does not exist");
    }

  //
  // Open file
  //
  this->m_InputGeometryFile = new std::ifstream;
  this->m_InputGeometryFile->open( this->m_FileName.c_str(), std::ios::in | std::ios::binary );

  if( this->m_InputGeometryFile->fail() )
    {
    this->ReleaseResources();
    itkExceptionMacro("Unable to open geometry file inputFilename= " << this->m_FileName );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::CloseGeometryFile()
{
  if( this->m_InputGeometryFile )
    {
    this->m_InputGeometryFile->close();
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::OpenDataFile()
{
  //
  // Open file
  //
  this->m_InputDataFile = new std::ifstream;
  this->m_InputDataFile->open( this->m_DataFileName.c_str(), std::ios::in | std::ios::binary );

  if( this->m_InputDataFile->fail() )
    {
    this->ReleaseResources();
    itkExceptionMacro("Unable to open data file inputFilename= " << this->m_DataFileName );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::CloseDataFile()
{
  this->m_InputDataFile->close();
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadGeometryHeader()
{
  this->ReadFileTypeFromGeometryFile();
  this->ReadCommentFromGeometryFile();
  this->ReadNumberOfPointsFromGeometryFile();
  this->ReadNumberOfCellsFromGeometryFile();
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadFileTypeFromGeometryFile()
{
  const unsigned int fileTypeIdLength = 3;
  char               fileTypeId[fileTypeIdLength];

  if( this->m_InputGeometryFile )
    {
    this->m_InputGeometryFile->read( fileTypeId, fileTypeIdLength );

    this->m_FileType = 0;
    this->m_FileType <<= 8;

    this->m_FileType |= fileTypeId[0];
    this->m_FileType <<= 8;

    this->m_FileType |= fileTypeId[1];
    this->m_FileType <<= 8;

    this->m_FileType |= fileTypeId[2];
    this->m_FileType <<= 8;

    itk::ByteSwapper<itk::uint32_t>::SwapFromSystemToBigEndian( &(this->m_FileType) );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadCommentFromGeometryFile()
{
  char byte;

  //
  //  Extract Comment, and ignore it.
  //
  byte = this->m_InputGeometryFile->get();

  this->m_Comment = "";

  unsigned int count = 0;

  while( byte != 0x0a && count < 45 )
    {
    this->m_Comment += byte;
    byte = this->m_InputGeometryFile->get();
    count++;
    }

  //
  // Try to get the second '\n', but if the '\n' is not there, we put the byte
  // back.
  //
  byte = this->m_InputGeometryFile->get();
  if( byte != 0x0a )
    {
    this->m_InputGeometryFile->unget();
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadNumberOfPointsFromGeometryFile()
{
  this->ReadInteger32( this->m_InputGeometryFile, this->m_NumberOfPoints );
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadNumberOfCellsFromGeometryFile()
{
  this->ReadInteger32( this->m_InputGeometryFile, this->m_NumberOfCells );
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadDataHeader()
{
  this->ReadFileTypeFromDataFile();
  this->ReadNumberOfPointsFromDataFile();
  this->ReadNumberOfCellsFromDataFile();
  this->ReadNumberOfValuesPerPointFromDataFile();
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadFileTypeFromDataFile()
{
  const unsigned int fileTypeIdLength = 3;
  char               fileTypeId[fileTypeIdLength];

  if( this->m_InputDataFile )
    {
    this->m_InputDataFile->read( fileTypeId, fileTypeIdLength );

    this->m_DataFileType = 0;
    this->m_DataFileType <<= 8;

    this->m_DataFileType |= fileTypeId[0];
    this->m_DataFileType <<= 8;

    this->m_DataFileType |= fileTypeId[1];
    this->m_DataFileType <<= 8;

    this->m_DataFileType |= fileTypeId[2];
    this->m_DataFileType <<= 8;

    itk::ByteSwapper<itk::uint32_t>::SwapFromSystemToBigEndian( &(this->m_DataFileType) );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadNumberOfPointsFromDataFile()
{
  itk::uint32_t numberOfPointsInDataFile;

  this->ReadInteger32( this->m_InputDataFile, numberOfPointsInDataFile );

  if( numberOfPointsInDataFile != this->m_NumberOfPoints )
    {
    this->ReleaseResources();
    itkExceptionMacro("Data file has " << numberOfPointsInDataFile
                                       << " points while Geometry file has " << this->m_NumberOfPoints );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadNumberOfCellsFromDataFile()
{
  itk::uint32_t numberOfCellsInDataFile;

  this->ReadInteger32( this->m_InputDataFile, numberOfCellsInDataFile );

  if( numberOfCellsInDataFile != this->m_NumberOfCells )
    {
    this->ReleaseResources();
    itkExceptionMacro("Data file has " << numberOfCellsInDataFile
                                       << " points while Geometry file has " << this->m_NumberOfCells );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadNumberOfValuesPerPointFromDataFile()
{
  itk::uint32_t numberOfValuesPerPoint;

  this->ReadInteger32( this->m_InputDataFile, numberOfValuesPerPoint );
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::PrepareOutputMesh()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  outputMesh->SetCellsAllocationMethod(
    OutputMeshType::CellsAllocatedDynamicallyCellByCell );
  outputMesh->GetPoints()->Reserve( this->m_NumberOfPoints );
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadSurface()
{
  this->ReadPoints();
  this->ReadCells();
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadPoints()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  PointType point;
  for( unsigned int ip = 0; ip < this->m_NumberOfPoints; ip++ )
    {
    this->ReadPoint( point );
    outputMesh->SetPoint( ip, point );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadCells()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  for( unsigned int ic = 0; ic < this->m_NumberOfCells; ic++ )
    {
    TriangleCellType * triangleCell = new TriangleCellType;

    this->ReadCell( *triangleCell );

    CellAutoPointer cell;
    cell.TakeOwnership( triangleCell );
    outputMesh->SetCell( ic, cell );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadPointData()
{
  typename PointDataContainer::Pointer pointData = PointDataContainer::New();
  pointData->Reserve( this->m_NumberOfPoints );

  PointDataIterator pointDataIterator = pointData->Begin();

  PointDataIterator endPoint = pointData->End();

  float pixelValue;

  while( pointDataIterator != endPoint )
    {
    this->m_InputDataFile->read( (char *)(&pixelValue), sizeof( pixelValue ) );
    itk::ByteSwapper<float>::SwapFromSystemToBigEndian( &pixelValue );
    pointDataIterator.Value() = static_cast<PixelType>( pixelValue );
    ++pointDataIterator;
    }

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  outputMesh->SetPointData( pointData );
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadCell( TriangleCellType & triangleCell )
{
  const unsigned int numberOfCellPoints = 3; // Triangles

  itk::uint32_t pointId;

  for( PointIdentifier k = 0; k < numberOfCellPoints; k++ )
    {
    this->ReadInteger32( this->m_InputGeometryFile, pointId );
    triangleCell.SetPointId( k, pointId );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadInteger32( std::ifstream * inputStream, itk::uint32_t & valueToRead )
{
  if( inputStream )
    {
    inputStream->read( (char *)(&valueToRead), sizeof(valueToRead) );
    itk::ByteSwapper<itk::uint32_t>::SwapFromSystemToBigEndian( &valueToRead );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadFloat( std::ifstream * inputStream, float & valueToRead  )
{
  if( inputStream )
    {
    inputStream->read( (char *)(&valueToRead), sizeof(valueToRead) );

    itk::ByteSwapper<float>::SwapFromSystemToBigEndian( &valueToRead );
    }
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::ReadPoint( PointType & point )
{
  float x;
  float y;
  float z;

  this->ReadFloat( this->m_InputGeometryFile, x );
  this->ReadFloat( this->m_InputGeometryFile, y );
  this->ReadFloat( this->m_InputGeometryFile, z );

  point[0] = x;
  point[1] = y;
  point[2] = z;
}

template <class TOutputMesh>
void
FreeSurferBinarySurfaceReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
  os << indent << "FileType: " << this->m_FileType << std::endl;
  os << indent << "Comment : " << this->m_Comment << std::endl;
  os << indent << "Number of Points : " << this->m_NumberOfPoints << std::endl;
  os << indent << "Number of Cells  : " << this->m_NumberOfCells << std::endl;

  os << indent << "DataFileName: " << this->m_DataFileName << std::endl;
  os << indent << "DataFileType: " << this->m_DataFileType << std::endl;
}
} // end of namespace itk

#endif
