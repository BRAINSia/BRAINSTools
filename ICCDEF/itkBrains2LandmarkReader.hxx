/*=================================================
  Program: itkBrains2LandmarkReader.h
  Date:    2009-10-22 10:26
  Version: 1.0
  Author:  Yongqiang Zhao
==================================================*/

// #ifndef __itkBrains2LandmarkReader_txx
// #define __itkBrains2LandmarkReader_txx

#include "itkBrains2LandmarkReader.h"
#include <fstream>

namespace itk
{
template <class TPixelType, unsigned Dimension>
Brains2LandmarkReader<TPixelType, Dimension>
::Brains2LandmarkReader()
{
  m_FileName = "";
  m_PointSet = InputPointSetType::New();
  m_PPointSet = PointSetType::New();
  m_ReferenceImage = ImageType::New();
}

template <class TPixelType, unsigned Dimension>
void Brains2LandmarkReader<TPixelType, Dimension>
::Update()
{
  if( m_FileName == "" )
    {
    itkExceptionMacro("No input FileName");
    return;
    }

  std::ifstream inputFile;
  inputFile.open(m_FileName.c_str() );

  if( inputFile.fail() )
    {
    itkExceptionMacro("Unable to open file\n"
                      "inputFilename= " << m_FileName );
    return;
    }

  std::string line;

  while( !inputFile.eof() )
    {
    std::getline(inputFile, line);
    if( line.find("LANDMARK_NUMBER_OF_POINTS:") != std::string::npos )
      {
      break;
      }
    }

  std::string pointLine(line, strlen("LANDMARK_NUMBER_OF_POINTS: "), line.length() );
  int         numberOfPoints = -1;

  if( sscanf(pointLine.c_str(), "%8d", &numberOfPoints) != 1 )
    {
    itkExceptionMacro("ERROR: Failed to read numberOfPoints\n"
                      "       pointLine= " << pointLine );
    return;
    }

  itkDebugMacro("numberOfPoints= " << numberOfPoints );

  if( numberOfPoints < 1 )
    {
    itkExceptionMacro("numberOfPoints < 1"
                      << "       numberOfPoints= " << numberOfPoints );
    return;
    }

  while( !inputFile.eof() && line.find("LANDMARK_X_SIZE:") == std::string::npos )
    {
    std::getline( inputFile, line );
    }

  itkDebugMacro( "X Dim:" << line );

  std::string xdimLine( line, strlen("LANDMARK_X_SIZE: "), line.length() );

  int xDim = 0;

  if( sscanf(xdimLine.c_str(), "%8d", &xDim) != 1 )
    {
    itkExceptionMacro("ERROR: Failed to read x Dimension\n"
                      "       xDim= " << xDim );
    return;
    }

  while( !inputFile.eof() && line.find("LANDMARK_X_RESOLUTION:") == std::string::npos )
    {
    std::getline( inputFile, line );
    }

  itkDebugMacro( "X Reso:" << line );

  std::string xresLine( line, strlen("LANDMARK_X_RESOLUTION: "), line.length() );

  float xRes = 0.0;

  if( sscanf(xresLine.c_str(), "%17f", &xRes) != 1 )
    {
    itkExceptionMacro("ERROR: Failed to read x resolution.\n"
                      "       xRes= " << xRes );
    return;
    }

  while( !inputFile.eof() && line.find("LANDMARK_Y_SIZE:") == std::string::npos )
    {
    std::getline( inputFile, line );
    }

  itkDebugMacro( "Y Dim:" << line );

  std::string ydimLine( line, strlen("LANDMARK_Y_SIZE: "), line.length() );

  int yDim = 0;

  if( sscanf(ydimLine.c_str(), "%8d", &yDim) != 1 )
    {
    itkExceptionMacro("ERROR: Failed to read y Dimension\n"
                      "       yDim= " << yDim );
    return;
    }

  while( !inputFile.eof() && line.find("LANDMARK_Y_RESOLUTION:") == std::string::npos )
    {
    std::getline( inputFile, line );
//    std::cout<<"xRes:"<<line<<std::endl;
    }

  itkDebugMacro( "Y Reso:" << line );

  std::string yresLine( line, strlen("LANDMARK_Y_RESOLUTION: "), line.length() );

  float yRes = 0.0;

  if( sscanf(yresLine.c_str(), "%17f", &yRes) != 1 )
    {
    itkExceptionMacro("ERROR: Failed to read y resolution.\n"
                      "       yRes= " << yRes );
    return;
    }

  while( !inputFile.eof() && line.find("LANDMARK_Z_SIZE:") == std::string::npos )
    {
    std::getline( inputFile, line );
    }

  itkDebugMacro( "Z Dim:" << line );

  std::string zdimLine( line, strlen("LANDMARK_Z_SIZE: "), line.length() );

  int zDim = 0;

  if( sscanf(zdimLine.c_str(), "%8d", &zDim) != 1 )
    {
    itkExceptionMacro("ERROR: Failed to read z Dimension\n"
                      "       zDim= " << zDim );
    return;
    }

  while( !inputFile.eof() && line.find("LANDMARK_Z_RESOLUTION:") == std::string::npos )
    {
    std::getline( inputFile, line );
    }

  itkDebugMacro( "Z Reso:" << line );

  float zRes = 0.0;

  if( sscanf(zdimLine.c_str(), "%17f", &zRes) != 1 )
    {
    itkExceptionMacro("ERROR: Failed to read z resolution.\n"
                      "       zRes= " << zRes );
    return;
    }

  while( !inputFile.eof() && line.find("IPL_HEADER_END") == std::string::npos )
    {
    std::getline(inputFile, line);
    // std::cout<<"end:"<<line<<std::endl;
    }

  typename InputPointSetType::PointType p;
  typename ImageType::PointType pp;
  typename ImageType::IndexType index;
  int   id;
  float x, y, z;
  for( int i = 0; i < numberOfPoints; i++ )
    {
    std::getline(inputFile, line);
    //   std::cout<<"2end:"<<line<<std::endl;
    sscanf(line.c_str(), "%8d %17f %17f %17f", &id, &x, &y, &z );
    //   std::cout<<id<<":  "<<x<<","<<y<<","<<z<<std::endl;
    if( x > xDim || x < 0 || y > yDim || y < 0 || z > zDim || z < 0 )
      {
      itkExceptionMacro("ERROR: Out of Image range.\n");
      return;
      }
    p[0] = (unsigned int)x;
    p[1] = (unsigned int)y;     // direction
    p[2] = (unsigned int)z;
    m_PointSet->SetPoint(i, p);
    index[0] = (unsigned int)x;
    index[1] = zDim - 1 - (unsigned int)z;
    index[2] = (unsigned int)y;
    m_ReferenceImage->TransformIndexToPhysicalPoint(index, pp);
    m_PPointSet->SetPoint(i, pp);
    std::cout << "Physical point:" << pp << std::endl;
    }
}
} // end of namespace itk
