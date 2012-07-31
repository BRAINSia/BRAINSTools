// Contributors: Eun Young (Regina) Kim and Joy Matsui
// Edit binary mask based on given landmark and direction.
//
#include "Slicer3LandmarkIO.h"

#include "BinaryMaskEditorBasedOnLandmarksCLP.h"

class PlaneLandmarksType
{
  // The definition depends on Slicer3LandmarkIO.
  // Point type is simply itk::Pointe<double, 3>
public:
  PointType A;
  PointType B;
  PointType C;
private:
  PointType normal;

  void SetNormal()
  {
    // Determine AB and AC vector components
    PointType AB;

    AB[0] = B[0] - A[0];
    AB[1] = B[1] - A[1];
    AB[2] = B[2] - A[2];

    PointType AC;
    AC[0] = C[0] - A[0];
    AC[1] = C[1] - A[1];
    AC[2] = C[2] - A[2];

    // Find cross product components

    normal[0] = ( AB[1] * AC[2] ) - ( AC[1] * AB[2] );
    normal[1] = -1 * ( ( AB[0] * AC[2] ) - ( AC[0] * AB[2] ) );
    normal[2] = ( AB[0] * AC[1] ) - ( AC[0] * AB[1] );
  }

  double GetRelativeLocationToPlane( PointType x )
  {
    double answer =
      normal[0] * ( A[0] - x[0] )
      + normal[1] * ( A[1] - x[1] )
      +  normal[2] * ( A[2] - x[2] );

    return answer;
  }
};

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // check input
  if( inputBinaryVolume.empty() ||
      inputLandmarksFilename.empty() ||
      inputLandmarkNames.empty() )
    {
    std::cout << "Input Landmarks Filename ( inputLandmarkFilename ) ,"
              << "input binary volume ( inputBinaryVolume ) and "
              << "input landmark name (inputLandmarkNames) are necessary."
              << std::endl;
    exit(EXIT_FAILURE);
    }
  if( inputLandmarkNames.size() != setCutDirectionForLandmark.size() )
    {
    std::cout << " Size should match between inputLandmarkNames and"
              << " setCutDirectionForLandmark but "
              << inputLandmarkNames.size() << " != "
              << setCutDirectionForLandmark.size() << std::endl;
    }

  // read inputBinaryVolume
  typedef unsigned char PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType, Dimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ImageReaderType;

  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( inputBinaryVolume );

  ImageType::Pointer inputVolume = imageReader->GetOutput();

  // read landmark file in
  std::cout << "Reading: "
            << inputLandmarksFilename << std::endl;
  LandmarksMapType landmarksSet = ReadSlicer3toITKLmk( inputLandmarksFilename );

  // duplicate image
  typdef itk::ImageDuplicator<ImageType> ImageDuplicatorType;
  ImageDuplicatorType::Pointer           duplicator = ImageDuplicatorType::New();
  duplicator->SetInputImage( inputVolume );
  duplicator->Update();

  ImageType::Pointer outputVolume =  duplicator->GetOutput();

  // cut by landmarks
  typedef std::vector<std::string>::const_iterator stringVectorIteratorType;
  for( stringVectorIteratorType ldmkIt = inputLandmarkNames.begin(),
       stringVectorIteratorType dircIt = setCutDirectionForLandmark.begin();
       ldmkIt < inputLandmarkNames.end(),
       ldmkIt++, dircIt++ )
    {
    if( landmarksSet.find( *ldmkIt ) == landmarksSet.end() )
      {
      std::cerr << "ERROR: Landmark not found: " << *ldmkIt << std::endl;
      std::exit( EXIT_FAILURE );
      }

    // TODO: implement CutBinaryVolumeByPointWithDirection
    outputVolume =
      CutBinaryVolumeByPointWithDirection( outputVolume, *ldmkIt, *dircIt );
    }

  // cut by plane defined by three landmarks

  typedef std::vector<std::string>           PlaneLandmarksType;
  typedef PlaneLandmarksType::const_iterator ThreeLandmarkIteratorType;

  typedef std::vector<PlaneLandmarksType>      PlaneLandmarkSetType;
  typedef PlaneLandmarkSetType::const_iterator PlaneLandmarkSetIteratorType;
  for( PlaneLandmarkSetIteratorType planeIt = inputLandmarkNamesForObliquePlane.begin();
       planeIt < inputLandmarkNamesForObliquePlane.end();
       planeIt++ )
    {
    PlaneLandmarksType currentPlane;

    ThreeLandmarkIteratorType tlmIt = ( *planeIt )->begin();
    for( unsigned int it = 0; it < 3; it++  ) // make sure we have only three ldmr per plane
      {
      if( landmarksSet.find( *tlmIt ) ==  landmarksSet.end() )
        {
        std::cerr << "ERROR: Landmark not found: " << *tlmIt << std::endl;
        std::exit( EXIT_FAILURE );
        }
      if( it == 1 )
        {
        currentPlane.A = *tlmIt;
        }
      if( it == 2 )
        {
        currentPlane.B = *tlmIt;
        }
      if( it == 3 )
        {
        currentPlane.C = *tlmIt;
        }

      tlmIt++;
      }
    currentPlane.SetNormal();

    outputVolume =
      CutBinaryVolumeByPlaneWithDirection( outputVolume, *currentPlane, *direction );
    }
}
