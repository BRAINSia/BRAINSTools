#include "Slicer3LandmarkIO.h"
#include "itkLandmarkBasedTransformInitializer.h"
#include "itkAffineTransform.h"
#include "itkImage.h"
#include "GenericTransformImage.h"

#include "LandmarkTransformInitializerCLP.h"

void CheckLandmarks( LandmarksMapType ldmk );

int
main(int argc, char * *argv)
{
  PARSE_ARGS;
  //  inputFixedLandmarkFilename
  //  inputMovingLandmarkFilename

  LandmarksMapType fixedLandmarks = ReadSlicer3toITKLmk( inputFixedLandmarkFilename );
  LandmarksMapType movingLandmarks = ReadSlicer3toITKLmk( inputMovingLandmarkFilename );

  typedef LandmarksMapType::const_iterator LandmarkIterator;
  for( LandmarkIterator it = fixedLandmarks.begin();
       it != fixedLandmarks.end();
       ++it )
    {
    std::cout << "fixed: " << it->first << " , " << it->second << std::endl;
    }

  /** check four landmarks */
  CheckLandmarks( fixedLandmarks );
  CheckLandmarks( movingLandmarks );

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
  for( LandmarkIterator fixedIt = fixedLandmarks.begin();
       fixedIt != fixedLandmarks.end();
       ++fixedIt )
    {
    LandmarkIterator movingIt = movingLandmarks.find( fixedIt->first );

    if( movingIt != movingLandmarks.end() )
      {
      fixedLmks.push_back( fixedIt->second);
      movingLmks.push_back( movingIt->second);
      }
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
