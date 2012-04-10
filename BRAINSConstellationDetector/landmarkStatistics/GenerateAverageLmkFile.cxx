/*
 * Author: Ali Ghayoor
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2011
 */

// This program gets n fcsv files as inputs. The first argument is the number of files. After that the fcsv file names
// come consequently.
// For instans if we have two landmark files as input:

//              .../GenerateAverageLmkFile 2 lmk1.fcsv lmk2.fcsv

// We have several fcsv file each one contains several landmarks with the same name but slightly different coordinates
// For EACH landmark we compute the average coordination.

// The final output of this program is a new landmark fcsv file which contains the average coordinate of each landmark

#include "itkImage.h"
#include "math.h"
#include "Slicer3LandmarkIO.h"

int main( int argc, char * argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: GenerateAverageLmkFile <number of files> <file1> <file2> ... <filen>"
              << std::endl;
    exit(EXIT_FAILURE);
    }

  const unsigned int k = atoi(argv[1]);  // The number of input landmark files

  if( argc != static_cast<int>(k + 2) )
    {
    std::cout << "First argument indicates the number of input landmarks files.\n"
              << "The number of input files inserted at the command line does not match." << std::endl;

    return EXIT_FAILURE;
    }

  unsigned int numNamedLandmarks = 0;
  double       x_ave, y_ave, z_ave;

  std::map<std::string, PointType> LandmarksAverageMap;

  typedef std::vector<std::map<std::string, PointType> > LandmarksMapTypeVec;
  LandmarksMapTypeVec LandmarksMapVector;

  // LandmarksMapType is as "std::map<std::string, PointType>" which means a map between landmarks and their
  // coordinates.
  // For each input landmark file this LandmarksMapType is computed and is set in a vector: "LandmarksMapTypeVec"

  LandmarksMapType temp;
  for( unsigned int i = 0; i < k; i++ )
    {
    temp = ReadSlicer3toITKLmk( argv[i + 2] );
    LandmarksMapVector.push_back(temp);
    }

  // Putting all the landmark names in a vector of strings and computing the number of landmarks
  std::vector<std::string> LandmarksNames;
  for( LandmarksMapType::const_iterator it = LandmarksMapVector[0].begin(); it != LandmarksMapVector[0].end(); ++it )
    {
    if( ( it->first ).compare("") != 0 )
      {
      LandmarksNames.push_back(it->first);
      numNamedLandmarks++;
      }
    }
  // Computing the average coordinate for each landmark
  for( unsigned int j = 0; j < numNamedLandmarks; j++ )
    {
    x_ave = y_ave = z_ave = 0;
    std::string name = LandmarksNames[j];
    for( unsigned int i = 0; i < k; i++ )
      {
      x_ave += LandmarksMapVector[i][name][0];
      y_ave += LandmarksMapVector[i][name][1];
      z_ave += LandmarksMapVector[i][name][2];
      }
    x_ave = x_ave / k;
    y_ave = y_ave / k;
    z_ave = z_ave / k;

    LandmarksAverageMap[name][0] = x_ave;
    LandmarksAverageMap[name][1] = y_ave;
    LandmarksAverageMap[name][2] = z_ave;
    }

  // writing the average of the landmarks in a new landmark file
  std::string initialFileName = argv[2];
  int         end = initialFileName.length() - 1;
  std::string outputLandmarkFileName = "a" + initialFileName.substr(1, end);
  // std::cout << outputLandmarkFileName << "is written to the disk." << std::endl;

  WriteITKtoSlicer3Lmk( outputLandmarkFileName, LandmarksAverageMap );

  return EXIT_SUCCESS;
}
