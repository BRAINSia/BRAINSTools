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
/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "Slicer3LandmarkIO.h"
#include "itkNumberToString.h"
#include "itkImageFileReaderException.h"

LandmarkWeightMapType ReadLandmarkWeights( const std::string & weightFilename )
{
  std::ifstream weightFileStream( weightFilename.c_str() );

  if( !weightFileStream.is_open() )
    {
    std::cerr << "Fail to open weight file " << std::endl;
    throw itk::ImageFileReaderException(__FILE__, __LINE__, "Couldn't open landmark weight file for reading", ITK_LOCATION);
    }

  std::string           line;
  LandmarkWeightMapType landmarkWeightMap;
  while( getline( weightFileStream, line ) )
    {
    const size_t      firstComma = line.find(',', 0);
    const std::string landmark = line.substr( 0, firstComma );
    const float       weight   = std::stod( (line.substr( firstComma + 1, line.length() - 1 ) ).c_str() );
    landmarkWeightMap[landmark] = weight;
    }

  return landmarkWeightMap;
}


static void WriteITKtoSlicer3LmkOldSlicer3( const std::string & landmarksFilename,
  const LandmarksMapType & landmarks )
{
  itk::NumberToString<double>    doubleConvert;
  const std::string fullPathLandmarksFileName = itksys::SystemTools::CollapseFullPath( landmarksFilename.c_str() );

  std::stringstream lmkPointStream;

  unsigned int numNamedLandmarks = 0;

  for( LandmarksMapType::const_iterator it = landmarks.begin(); it != landmarks.end(); ++it )
  {
    if( ( it->first ).compare("") != 0 )
    {
      // NOTE: Slicer3 use RAS coordinate system to represent landmarks
      // but ITK landmarks are in LPS, so we need to negate the first two
      // component of the landmark points.
      lmkPointStream <<     it->first       << ","
                     << doubleConvert(-( it->second[0] ) ) << ","
                     << doubleConvert(-( it->second[1] ) ) << ","
                     << doubleConvert(+( it->second[2] ) )
                     << ",1,1\n"; // Note the last two columns are
      // ,visible,editable
      ++numNamedLandmarks;
    }
  }

  std::stringstream lmksStream;
  // Write the .fcvs header information.
  lmksStream << "#Fiducial List file " << fullPathLandmarksFileName << std::endl;
  lmksStream << "#numPoints = " << numNamedLandmarks << "\n";
  lmksStream << "#symbolScale = 5" << std::endl;
  lmksStream << "#visibility = 1" << std::endl;
  lmksStream << "#textScale = 4.5" << std::endl;
  lmksStream << "#color = 0.4,1,1" << std::endl;
  lmksStream << "#selectedColor = 1,0.5,0.5" << std::endl;
  lmksStream << "#label,x,y,z,sel,vis" << std::endl;
  lmksStream << lmkPointStream.str();

  // Now write file to disk
  std::ofstream myfile;
  myfile.open( fullPathLandmarksFileName.c_str() );
  if( !myfile.is_open() )
  {
    std::cerr << "Error: Can't write Slicer3 landmark file "
              << fullPathLandmarksFileName << std::endl;
    std::cerr.flush();
    throw itk::ImageFileReaderException(__FILE__, __LINE__, "Couldn't open landmarks file for writing", ITK_LOCATION);
  }
  myfile << lmksStream.str();
  myfile.close();
}


static void WriteITKtoSlicer3LmkSlicer4( const std::string & landmarksFilename,
  const LandmarksMapType & landmarks )
{
  itk::NumberToString<double>    doubleConvert;
  const std::string fullPathLandmarksFileName = itksys::SystemTools::CollapseFullPath( landmarksFilename.c_str() );

  std::stringstream lmkPointStream;
  unsigned int numNamedLandmarks = 0;
  for( LandmarksMapType::const_iterator it = landmarks.begin(); it != landmarks.end(); ++it )
  {
    if( ( it->first ).compare("") != 0 )
    {
      // NOTE: Slicer4 Markups use RAS coordinate system to represent landmarks
      // but ITK landmarks are in LPS, so we need to negate the first two
      // component of the landmark points.
      lmkPointStream << "vtkMRMLMarkupsFiducialNode_" <<numNamedLandmarks << ",";

      lmkPointStream << doubleConvert(-( it->second[0] ) ) << ","
                     << doubleConvert(-( it->second[1] ) ) << ","
                     << doubleConvert(+( it->second[2] ) )
                     << ",0,0,0,1,1,1,0,"<< it->first << ",," << std::endl; // Note the last two columns are
      ++numNamedLandmarks;
    }
  }

  std::stringstream lmksStream;
  // Write the .fcvs header information.
  lmksStream << "# Markups fiducial file version = 4.6" << std::endl;
  lmksStream << "# CoordinateSystem = 0" << std::endl;
  lmksStream << "# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID" << std::endl;
  lmksStream << lmkPointStream.str();

  // Now write file to disk
  std::ofstream myfile;
  myfile.open( fullPathLandmarksFileName.c_str() );
  if( !myfile.is_open() )
  {
    std::cerr << "Error: Can't write Slicer3 landmark file "
              << fullPathLandmarksFileName << std::endl;
    std::cerr.flush();
    throw itk::ImageFileReaderException(__FILE__, __LINE__, "Couldn't open landmarks file for writing", ITK_LOCATION);
  }
  myfile << lmksStream.str();
  myfile.close();
}


static LandmarksMapType
ReadSlicer3toITKLmkOldSlicer( const std::string & landmarksFilename )
{
  LandmarksMapType landmarks;
  std::string      landmarksFilenameTmp = itksys::SystemTools::CollapseFullPath( landmarksFilename.c_str() );
  std::ifstream    myfile( landmarksFilenameTmp.c_str() );

  if( !myfile.is_open() )
  {
    std::cerr << "Error: Failed to load landmarks file!" << std::endl;
    std::cerr.flush();
    std::string errorMsg = std::string("Couldn't open landmarks file for reading: ")+landmarksFilename;
    throw itk::ImageFileReaderException(__FILE__, __LINE__, errorMsg.c_str(), ITK_LOCATION);
    // do not return empty landmarks
  }
  std::string line;
  while( getline( myfile, line ) )
  {
    if( line.compare( 0, 1, "#" ) != 0 )  // Skip lines starting with a #
    {
      size_t            pos1 = line.find( ',', 0 );
      const std::string name = line.substr( 0, pos1 );
      LandmarkPointType labelPos;
      for( unsigned int i = 0; i < 3; ++i )
      {
        const size_t pos2 = line.find( ',', pos1 + 1 );
        labelPos[i] = std::stod( line.substr(pos1 + 1, pos2 - pos1 - 1 ).c_str() );
        if( i < 2 )  // Negate first two components for RAS -> LPS
        {
          labelPos[i] *= -1;
        }
        pos1 = pos2;
      }
      landmarks[name] = labelPos;
    }
  }

  myfile.close();
  return landmarks;
}

static std::vector<std::string> split(const std::string& s, char delim)
{
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> tokens;
  while (getline(ss, item, delim)) {
    tokens.push_back(item);
  }
  return tokens;
}

static LandmarksMapType
ReadSlicer3toITKLmkSlicer4(const std::string& landmarksFilename)
{
  LandmarksMapType landmarks;

  std::string landmarksFilenameTmp = itksys::SystemTools::CollapseFullPath(landmarksFilename.c_str());
  std::ifstream myfile(landmarksFilenameTmp.c_str());

  if (!myfile.is_open()) {
    std::cerr << "Error: Failed to load landmarks file!" << std::endl;
    std::cerr.flush();
    std::string errorMsg = std::string("Couldn't open landmarks file for reading: ")+landmarksFilename;
    throw itk::ImageFileReaderException(__FILE__, __LINE__, errorMsg.c_str(), ITK_LOCATION);
    // do not return empty landmarks
  }
  std::string line;
  while (getline(myfile, line)) {
    if (line.compare(0, 1, "#")!=0)  // Skip lines starting with a #
    {
      //             0 1 2 3  4  5  6  7   8   9   10    11   12               13
      //# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID
      std::vector<std::string> asTokens = split(line, ',');
      const std::string name = asTokens[11]; //11 position is label name
      LandmarkPointType labelPos;
      labelPos[0] = -atof(asTokens[1].c_str()); // L -> R (-)
      labelPos[1] = -atof(asTokens[2].c_str()); // P -> A (-)
      labelPos[2] = +atof(asTokens[3].c_str()); // S -> S (+)
      landmarks[name] = labelPos;
    }
  }
  myfile.close();
  return landmarks;
}


LandmarksMapType
ReadSlicer3toITKLmk( const std::string & landmarksFilename )
{
  std::ifstream    myfile( landmarksFilename.c_str() );
  if( !myfile.is_open() )
  {
    std::cerr << "Error: Failed to load landmarks file!" << std::endl;
    std::cerr.flush();
    std::string errorMsg = std::string("Couldn't open landmarks file for reading: ")+landmarksFilename;
    throw itk::ImageFileReaderException(__FILE__, __LINE__, errorMsg.c_str(), ITK_LOCATION);
    // do not return empty landmarks
  }
  std::string firstLine;
  std::getline(myfile,firstLine);
  myfile.close();
  if( firstLine.find("Fiducial List") != std::string::npos)
  {
    const LandmarksMapType landmarks =  ReadSlicer3toITKLmkOldSlicer(landmarksFilename);
    return landmarks;
  }
  else if(firstLine.find("Markups fiducial file") != std::string::npos )
  {
    const LandmarksMapType landmarks =  ReadSlicer3toITKLmkSlicer4(landmarksFilename);
    return landmarks;
  }
  else
  {
    itkGenericExceptionMacro(<< "ERROR:  Unrecognized landmark file format for " << landmarksFilename);
  }
}

void WriteITKtoSlicer3Lmk(const std::string& landmarksFilename,
  const LandmarksMapType& landmarks, const SLICER_LANDMARK_FILE_TYPE slicerLmkType)
{
  switch (slicerLmkType) {
  case SLICER_V3_FCSV:
    WriteITKtoSlicer3LmkOldSlicer3(landmarksFilename, landmarks);
    break;
  case SLICER_V4_FCSV:
    WriteITKtoSlicer3LmkSlicer4(landmarksFilename, landmarks);
    break;
  case SLICER_FCSV_BEGIN: itkGenericExceptionMacro(
    << "ERROR:  MODE NOT IMPLEMENTED FOR LANDMARK WRITING (SLICER_FCSV_BEGIN)");
    ITK_FALLTHROUGH;
  case SLICER_FCSV_END: itkGenericExceptionMacro(
    << "ERROR:  MODE NOT IMPLEMENTED FOR LANDMARK WRITING (SLICER_FCSV_END)");
    ITK_FALLTHROUGH;
  default:
    /* Do Nothing */
  itkGenericExceptionMacro(<< "ERROR:  MODE NOT IMPLEMENTED FOR LANDMARK WRITING (No Match): " << slicerLmkType);
    break;
  }
  return;
}
