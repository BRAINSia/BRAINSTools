/*
 * Author: Ali Ghayoor
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2013
 */

/*
 This program program uses the input transform to propagate
 reference landmarks to the target landmark file.
*/

#include "itkPoint.h"
#include "itkTransformFileReader.h"
#include "itkCompositeTransform.h"
#include "Slicer3LandmarkIO.h"

#include "BRAINSConstellationLandmarksTransformCLP.h"


int main( int argc, char *argv[] )
{
  PARSE_ARGS;

  if( inputLandmarksFile == ""
     || inputTransformFile == ""
     || outputLandmarksFile == "")
    {
    std::cerr << "Input and output file names should be given by commandline. " << std::endl;
    std::cerr << "Usage:\n"
              << "~/BRAINSConstellationLandmarksTransform\n"
              << "--inputLandmarksFile (-i)\n"
              << "--inputTransformFile (-t)\n"
              << "--outputLandmarksFile (-o)"
              << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;

  typedef itk::Point<double,Dimension> PointType;
  typedef std::map< std::string, PointType > LandmarksMapType;

  LandmarksMapType origLandmarks = ReadSlicer3toITKLmk( inputLandmarksFile );
  LandmarksMapType transformedLandmarks;

  typedef itk::TransformFileReader ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputTransformFile );
  reader->Update();

  ReaderType::TransformListType *transformList = reader->GetTransformList();

  typedef itk::CompositeTransform<double, Dimension> CompositeTransformType;

  CompositeTransformType::Pointer inputCompTrans =
  dynamic_cast<CompositeTransformType *>( transformList->front().GetPointer() );
  if( inputCompTrans.IsNull() )
    {
    std::cerr << "The input transform should be a composite transform." << std::endl;
    return EXIT_FAILURE;
    }

  LandmarksMapType::const_iterator it = origLandmarks.begin();
  for(; it!=origLandmarks.end(); it++)
    {
    transformedLandmarks[it->first] = inputCompTrans->TransformPoint( it->second );
    }

  WriteITKtoSlicer3Lmk( outputLandmarksFile, transformedLandmarks );
  std::cout << "The transformed landmarks file is written." << std::endl;

  return EXIT_SUCCESS;
}
