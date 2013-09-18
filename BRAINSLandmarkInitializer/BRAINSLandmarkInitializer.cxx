#include "itkLandmarkBasedTransformInitializer.h"
#include "itkAffineTransform.h"
#include "itkImage.h"

#include <map>

#include "GenericTransformImage.h"
#include "Slicer3LandmarkIO.h"

#include "BRAINSLandmarkInitializerCLP.h"


static void CheckLandmarks( const LandmarksMapType & ldmk, const LandmarkWeightMapType & weightMap)
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
  for( std::map<std::string, float>::const_iterator i = weightMap.begin();
       i != weightMap.end();
       ++i )
    {
    if( ldmk.find( i->first ) == ldmk.end() )
      {
      std::cerr << "WARNING: Landmark not found: " << i->first << std::endl;
      }
#if defined(VERBOSE_OUTPUT)
    else
      {
      std::cerr << "NOTE: Landmark found: " << i->first << std::endl;
      }
#endif
    }
}


int
main(int argc, char *argv[])
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

  /** Landmark Weights */
  LandmarkWeightMapType landmarkWeightMap = ReadLandmarkWeights( inputWeightFilename );

  /** read in *fcsv file */
  /** check four landmarks */
  std::cout << "Reading: " << inputFixedLandmarkFilename << std::endl;
  LandmarksMapType fixedLandmarks = ReadSlicer3toITKLmk( inputFixedLandmarkFilename );
  CheckLandmarks( fixedLandmarks, landmarkWeightMap );

  std::cout << "Reading: " << inputMovingLandmarkFilename << std::endl;
  LandmarksMapType movingLandmarks = ReadSlicer3toITKLmk( inputMovingLandmarkFilename );
  CheckLandmarks( movingLandmarks, landmarkWeightMap );

  /** Landmark Initializaer */
  typedef double PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType, Dimension>           ImageType;
  typedef itk::AffineTransform<PixelType, Dimension> AffineTransformType;
  AffineTransformType::Pointer affineTransform = AffineTransformType::New();

  typedef itk::LandmarkBasedTransformInitializer<AffineTransformType,
                                                 ImageType,
                                                 ImageType> LandmarkBasedInitializerType;

  LandmarkBasedInitializerType::Pointer landmarkBasedInitializer =
    LandmarkBasedInitializerType::New();

  typedef LandmarkBasedInitializerType::LandmarkWeightType LandmarkWeightContainerType;
  LandmarkWeightContainerType landmarkWgts;

  typedef LandmarkBasedInitializerType::LandmarkPointContainer LandmarkContainerType;
  LandmarkContainerType fixedLmks;
  LandmarkContainerType movingLmks;
  typedef LandmarksMapType::const_iterator LandmarkConstIterator;
  for( LandmarkConstIterator fixedIt = fixedLandmarks.begin();
       fixedIt != fixedLandmarks.end();
       ++fixedIt )
    {
    LandmarkConstIterator movingIt = movingLandmarks.find( fixedIt->first );
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

  itk::WriteTransformToDisk<double>( affineTransform, outputTransformFilename);

  return EXIT_SUCCESS;
}
