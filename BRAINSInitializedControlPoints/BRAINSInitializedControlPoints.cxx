/*=========================================================================

Module:    $RCSfile: BRAINSInitializedControlPoints.cxx,v $
Language:  C++
Date:      $Date: 2010/11/29 10:37:07 $
Version:   $Revision: 0.1 $

This work is part of the National Alliance for Medical Image
Computing (NAMIC), funded by the National Institutes of Health
through the NIH Roadmap for Medical Research, Grant U54 EB005149.

See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

 ***
 This program converts takes an image and the bspline grid spacing then
 outputs the bspline control points as a landmark file.

 =========================================================================*/

#include <sstream>
#include "itkXMLFilterWatcher.h"

#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itksys/Base64.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorImage.h"

#include "itkBSplineDeformableTransform.h"
#include "itkBSplineDeformableTransformInitializer.h"
#include "Slicer3LandmarkIO.h"

#include "BRAINSInitializedControlPointsCLP.h"
#include "BRAINSThreadControl.h"

#include "itkPermuteAxesImageFilter.h"

int main(int argc, char* argv[])
{
  PARSE_ARGS;
  BRAINSUtils::SetThreadCount(numberOfThreads);

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
    }
  if( splineGridSize.size() == 0 )
    {
    violated = true; std::cout << "  --splineGridSize Required! "  << std::endl;
    }
  if( outputLandmarksFile.size() == 0 )
    {
    violated = true; std::cout << "  --outputLandMarksFile Required! "  << std::endl;
    }
  if( permuteOrder.size() != 3 )
    {
    violated = true; std::cout << " --permuteOrder must have size 3! " << std::endl;
    }
  if( violated )
    {
    exit(1);
    }

  typedef float PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType, Dimension> FixedVolumeType;

  // Verify that the spline grid sizes are greater than 3
    {
    for( unsigned int sgs = 0; sgs < splineGridSize.size(); sgs++ )
      {
      if( splineGridSize[sgs] < 3 )
        {
        std::cout << "splineGridSize[" << sgs << "]= " << splineGridSize[sgs]
                  << " is invalid.  There must be at lest 3 divisions in each dimension of the image." << std::endl;
        exit(-1);
        }
      }
    }

  FixedVolumeType::Pointer fixedImage = NULL;
  try
    {
    typedef itk::ImageFileReader<FixedVolumeType> FixedVolumeReaderType;
    FixedVolumeReaderType::Pointer fixedVolumeReader = FixedVolumeReaderType::New();
    fixedVolumeReader->SetFileName(inputVolume);
    fixedVolumeReader->Update();
    typedef itk::PermuteAxesImageFilter<FixedVolumeType> PermuterType;
    PermuterType::PermuteOrderArrayType myOrdering;
    myOrdering[0] = permuteOrder[0];
    myOrdering[1] = permuteOrder[1];
    myOrdering[2] = permuteOrder[2];
    PermuterType::Pointer myPermuter = PermuterType::New();
    myPermuter->SetInput(fixedVolumeReader->GetOutput() );
    myPermuter->SetOrder(myOrdering);
    myPermuter->Update();
    fixedImage = myPermuter->GetOutput();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::BSplineDeformableTransform<float, 3, 3> BSplineTransformType;
  BSplineTransformType::Pointer initialBSplineTransform = BSplineTransformType::New();
  initialBSplineTransform->SetIdentity();

  typedef BSplineTransformType::RegionType TransformRegionType;
  typedef TransformRegionType::SizeType    TransformSizeType;

  typedef itk::BSplineDeformableTransformInitializer
    <BSplineTransformType, FixedVolumeType> InitializerType;
  InitializerType::Pointer transformInitializer = InitializerType::New();

  transformInitializer->SetTransform(initialBSplineTransform);
  transformInitializer->SetImage(fixedImage);
  TransformSizeType tempGridSize;
  tempGridSize[0] = splineGridSize[0];
  tempGridSize[1] = splineGridSize[1];
  tempGridSize[2] = splineGridSize[2];

  std::cout << "splineGridSize=" << splineGridSize[0] << "," << splineGridSize[1] << "," << splineGridSize[2]
            << std::endl;

  transformInitializer->SetGridSizeInsideTheImage(tempGridSize);
  try
    {
    transformInitializer->InitializeTransform();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  // HACK:   The following silly setting is necessary due to a race condition
  // that
  //        is present in the BSplineTransformationCode between the differing
  //        representations of the parameters
  BSplineTransformType::ParametersType temp = initialBSplineTransform->GetParameters();
  initialBSplineTransform->SetParameters(temp);
  // std::cout << transformInitializer << std::endl;

  BSplineTransformType::ImagePointer coefficientImage = initialBSplineTransform->GetCoefficientImage()[0];
  std::cout << "fixedImage\n" << fixedImage << std::endl;
  std::cout << "coefficientImage\n" << coefficientImage << std::endl;
  std::cout << "initialBSplineTransform\n" << initialBSplineTransform << std::endl;

  typedef itk::ImageRegionIterator<BSplineTransformType::ImageType> IteratorType;

  IteratorType coefItr( coefficientImage, coefficientImage->GetLargestPossibleRegion() );
  coefItr.GoToBegin();

  LandmarksMapType outputLandmarks;
  unsigned int     locationCount = 0;
  while( !coefItr.IsAtEnd() )
    {
    BSplineTransformType::ImageType::PointType currentPoint;
    coefficientImage->TransformIndexToPhysicalPoint(coefItr.GetIndex(), currentPoint);
    std::stringstream tmpName("p");
    tmpName << locationCount;
    locationCount++;
    outputLandmarks[tmpName.str()] = currentPoint;
    ++coefItr;
    }

  typedef itk::ImageFileWriter<FixedVolumeType> FixedVolumeWriterType;
  FixedVolumeWriterType::Pointer permutedWriter = FixedVolumeWriterType::New();
  permutedWriter->SetInput(fixedImage);
  permutedWriter->SetFileName(outputVolume);
  permutedWriter->Update();

  WriteITKtoSlicer3Lmk( outputLandmarksFile, outputLandmarks);
}
