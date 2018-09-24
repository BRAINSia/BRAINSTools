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
 * Author: Ali Ghayoor, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2011
 */

#include "Slicer3LandmarkWeightIO.h"
#include "itkNumberToString.h"
void
WriteITKtoSlicer3LmkWts( const std::string & landmarksWeightFilename,
                         const LandmarksWeightMapType & landmarks )
{
  itk::NumberToString<double> doubleConvert;

  const std::string fullPathLandmarksWeightFileName =
    itksys::SystemTools::CollapseFullPath( landmarksWeightFilename.c_str() );

  std::stringstream lmkWeightStream;

  unsigned int numNamedLandmarks = 0;

  for( LandmarksWeightMapType::const_iterator it = landmarks.begin(); it != landmarks.end(); ++it )
    {
    if( ( it->first ).compare("") != 0 )
      {
      lmkWeightStream << it->first << ","
                      << doubleConvert(it->second) << std::endl;

      ++numNamedLandmarks;
      }
    }

  std::stringstream lmksStream;
  // Write the .wts header information.
//  lmksStream << "#Fiducial List file " << fullPathLandmarksWeightFileName << std::endl;
//  lmksStream << "#numPoints = " << numNamedLandmarks << "\n";
//  lmksStream << "#symbolScale = 5" << std::endl;
//  lmksStream << "#visibility = 1" << std::endl;
//  lmksStream << "#textScale = 4.5" << std::endl;
//  lmksStream << "#color = 0.4,1,1" << std::endl;
//  lmksStream << "#selectedColor = 1,0.5,0.5" << std::endl;
  lmksStream << "#label,weight" << std::endl;
  lmksStream << lmkWeightStream.str();

  // Now write file to disk
  std::ofstream myfile;
  myfile.open( fullPathLandmarksWeightFileName.c_str() );
  if( !myfile.is_open() )
    {
    std::cerr << "Error: Can't write Slicer3 landmark weight file "
              << fullPathLandmarksWeightFileName << std::endl;
    std::cerr.flush();
    return;
    }
  myfile << lmksStream.str();
  myfile.close();
}

extern LandmarksWeightMapType
ReadSlicer3toITKLmkWts( const std::string & landmarksWeightFilename )
{
  LandmarksWeightMapType landmarks;
  std::string            landmarksFilenameTmp =
    itksys::SystemTools::CollapseFullPath( landmarksWeightFilename.c_str() );
  std::ifstream myfile( landmarksFilenameTmp.c_str() );

  if( !myfile.is_open() )
    {
    std::cerr << "Error: Failed to load landmark weights file!" << std::endl;
    std::cerr.flush();
    return landmarks; // return empty landmarks
    }
  std::string line;
  while( getline( myfile, line ) )
    {
    if( line.compare( 0, 1, "#" ) != 0 )  // Skip lines starting with a #
      {
      size_t            pos1 = line.find( ',', 0 );
      const std::string name = line.substr( 0, pos1 );
      double            weight;
      const size_t      pos2 = line.find( ' ', pos1 + 1 );
      weight = std::stod( line.substr(pos1 + 1, pos2 - pos1 - 1 ).c_str() );

      landmarks[name] = weight;
      }
    }

  myfile.close();
  return landmarks;
}
