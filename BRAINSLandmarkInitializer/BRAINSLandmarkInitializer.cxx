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
#include "itkLandmarkBasedTransformInitializer.h"
#include "itkVersorRigid3DTransform.h"
#include "itkAffineTransform.h"
#include "itkBSplineTransform.h"
#include "itkImage.h"
#include "itkImageFileReader.h"

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

template<class TTransformType>
int
InitializeTransform( int argc, char *argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

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

  ImageType::Pointer referenceImage = nullptr;
  if( !inputReferenceImageFilename.empty() )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( inputReferenceImageFilename );
    try
      {
      reader->Update();
      referenceImage = reader->GetOutput();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cout << "Exception Object caught: " << std::endl;
      std::cout << err << std::endl;
      throw;
      }
    }

  typedef TTransformType                             LocalTransformType;
  typename LocalTransformType::Pointer transform =
    LocalTransformType::New();

  typedef itk::LandmarkBasedTransformInitializer<LocalTransformType,
                                                 ImageType,
                                                 ImageType> LandmarkBasedInitializerType;

  typename LandmarkBasedInitializerType::Pointer landmarkBasedInitializer =
    LandmarkBasedInitializerType::New();

  typedef typename LandmarkBasedInitializerType::LandmarkWeightType LandmarkWeightContainerType;
  LandmarkWeightContainerType landmarkWgts;

  typedef typename LandmarkBasedInitializerType::LandmarkPointContainer LandmarkContainerType;
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

  if( referenceImage.IsNotNull() )
    {
    landmarkBasedInitializer->SetReferenceImage( referenceImage );
    }
  landmarkBasedInitializer->SetBSplineNumberOfControlPoints( bsplineNumberOfControlPoints );

  landmarkBasedInitializer->SetFixedLandmarks( fixedLmks );
  landmarkBasedInitializer->SetMovingLandmarks( movingLmks);
  landmarkBasedInitializer->SetTransform( transform );
  try
    {
    landmarkBasedInitializer->InitializeTransform();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Exception Object caught: " << std::endl;
    std::cout << err << std::endl;
    throw;
    }

  itk::WriteTransformToDisk<double>( transform, outputTransformFilename);

  return EXIT_SUCCESS;
}

//////////////////// M A I N /////////////////////////////////////////////////
int
main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

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

  typedef double ParameterValueType;
  const unsigned int Dimension = 3;

  if( outputTransformType == "AffineTransform" )
    {
    typedef itk::AffineTransform<ParameterValueType, Dimension> AffineTransformType;
    return InitializeTransform<AffineTransformType>( argc, argv );
    }
  else if( outputTransformType == "BSplineTransform" )
    {
    const static unsigned int SplineOrder = 3;
    typedef itk::BSplineTransform< ParameterValueType,
                                   Dimension,
                                   SplineOrder>                  BSplineTransformType;
    return InitializeTransform<BSplineTransformType>( argc, argv );
    }
  else if( outputTransformType == "VersorRigid3DTransform" )
    {
    typedef itk::VersorRigid3DTransform< ParameterValueType >    VersorRigid3DTransformType;
    return InitializeTransform<VersorRigid3DTransformType>( argc, argv );
    }
  else
    {
    std::cerr << "Error: Invalid parameter for output transform type." << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
