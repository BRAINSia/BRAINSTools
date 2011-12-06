#include "itkLandmarkBasedTransformInitializer.h"
#include "itkAffineTransform.h"
#include "itkImage.h"

#include <map>

#include "GenericTransformImage.h"
#include "Slicer3LandmarkIO.h"

#include "BRAINSLandmarkInitializerCLP.h"

typedef std::map<std::string, float> LandmarkWeightMapType;

void CheckLandmarks( LandmarksMapType ldmk );

LandmarkWeightMapType ReadLandmarkWeights( std::string weightFilename );

int
main(int argc, char * *argv)
{
  PARSE_ARGS;

  if( inputFixedLandmarkFilename.empty() ||
      inputMovingLandmarkFilename.empty() ||
      outputTransformFilename.empty() )
    {
    std::cout << "Input Landmarks ( inputFixedLandmarkFilename ,"
              << "inputMovingLandmarkFilename ) and "
              << "outputTransformationFilename are necessary"
              << std::endl;
    exit(EXIT_FAILURE);
    }

  /** read in *fcsv file */
  LandmarksMapType fixedLandmarks = ReadSlicer3toITKLmk( inputFixedLandmarkFilename );
  LandmarksMapType movingLandmarks = ReadSlicer3toITKLmk( inputMovingLandmarkFilename );

  /** check four landmarks */
  CheckLandmarks( fixedLandmarks );
  CheckLandmarks( movingLandmarks );

  /** Landmark Weights */
  LandmarkWeightMapType landmarkWeightMap = ReadLandmarkWeights( inputWeightFilename );

  /** Landmark Initializaer */
  typedef double PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType, Dimension> ImageType;

  typedef itk::AffineTransform<PixelType, Dimension> AffineTransformType;
  AffineTransformType::Pointer affineTransform = AffineTransformType::New();

  typedef itk::LandmarkBasedTransformInitializer<AffineTransformType,
                                                 ImageType,
                                                 ImageType> LandmarkBasedInitializerType;

  LandmarkBasedInitializerType::Pointer landmarkBasedInitializer =
    LandmarkBasedInitializerType::New();

  typedef LandmarkBasedInitializerType::LandmarkPointContainer LandmarkContainerType;

  LandmarkContainerType fixedLmks;
  LandmarkContainerType movingLmks;

  typedef LandmarkBasedInitializerType::LandmarkWeightType LandmarkWeightContainerType;
  LandmarkWeightContainerType landmarkWgts;

  typedef LandmarksMapType::const_iterator LandmarkIterator;
  for( LandmarkIterator fixedIt = fixedLandmarks.begin();
       fixedIt != fixedLandmarks.end();
       ++fixedIt )
    {
    LandmarkIterator movingIt = movingLandmarks.find( fixedIt->first );

    if( movingIt != movingLandmarks.end() )
      {
      fixedLmks.push_back( fixedIt->second);
      movingLmks.push_back( movingIt->second);

      if( !inputWeightFilename.empty() )
        {
        if( landmarkWeightMap.find( fixedIt->first ) != landmarkWeightMap.end() )
          {
          landmarkWgts.push_back( landmarkWeightMap[fixedIt->first] );
          }
        else
          {
          std::cout << "Landmark for " << fixedIt->first << " does not exist. "
                    << "Set the weight to 0.5 "
                    << std::endl;
          landmarkWgts.push_back( 0.5F );
          }
        }
      }
    }
  /** set weights */
  if( !inputWeightFilename.empty() )
    {
    landmarkBasedInitializer->SetLandmarkWeight( landmarkWgts );
    }

  landmarkBasedInitializer->SetFixedLandmarks( fixedLmks );
  landmarkBasedInitializer->SetMovingLandmarks( movingLmks);

  landmarkBasedInitializer->SetTransform( affineTransform );

  landmarkBasedInitializer->InitializeTransform();

  WriteTransformToDisk( affineTransform,
                        outputTransformFilename);

  return EXIT_SUCCESS;
}

void
CheckLandmarks( LandmarksMapType ldmk )
{
  if( ldmk.size() < 4 )
    {
    std::cerr << "At least 3 fiducual points must be specified. "
              << std::endl;
    exit(EXIT_FAILURE);
    }

  if( ldmk.find( "AC" ) == ldmk.end() ||
      ldmk.find( "PC" ) == ldmk.end() ||
      ldmk.find( "LE" ) == ldmk.end() ||
      ldmk.find( "RE" ) == ldmk.end() )
    {
    std::cerr << " Base four landmarks ( AC, PC, left eye(LE), and right eye(RE) ) "
              << " has to be provided"
              << std::endl;
    exit(EXIT_FAILURE);
    }
}

LandmarkWeightMapType
ReadLandmarkWeights( std::string weightFilename )
{
  std::ifstream weightFileStream( weightFilename.c_str() );

  if( !weightFileStream.is_open() )
    {
    std::cerr << "Fail to open weight file "
              << std::endl;
    exit(EXIT_FAILURE);
    }

  std::string           line;
  LandmarkWeightMapType landmarkWeightMap;

  while( getline( weightFileStream, line ) )
    {
    size_t            firstComma = line.find(',', 0);
    const std::string landmark = line.substr( 0, firstComma );
    const float       weight   = atof( (line.substr( firstComma + 1, line.length() - 1 ) ).c_str() );
    landmarkWeightMap[landmark] = weight;
    }

  return landmarkWeightMap;
}
