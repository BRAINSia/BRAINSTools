/*
 * Author: Ali Ghayoor, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2011
 */

#include "Slicer3LandmarkWeightIO.h"
#include "DoubleToString.h"
void
WriteITKtoSlicer3LmkWts( const std::string & landmarksWeightFilename,
                         const LandmarksWeightMapType & landmarks )
{
  DoubleToString doubleConvert;

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
      weight = atof( line.substr(pos1 + 1, pos2 - pos1 - 1 ).c_str() );

      landmarks[name] = weight;
      }
    }

  myfile.close();
  return landmarks;
}
