/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab, University of Iowa Health Care, 2010
 */

#include "itkImage.h"
#include "itkPoint.h"
#include "itkPointSet.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkTransformFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "BRAINSThreadControl.h"

// Use modified itkKernelTransform to get affine transform
#include "itkKernelTransform.h"
#include "itkThinPlateSplineKernelTransform.h"

template <class TScalarType,
          unsigned int NDimension>
class BCDThinPlateSplineKernelTransform :
  public itk::ThinPlateSplineKernelTransform<TScalarType, NDimension>
{
public:
  /** Standard class typedefs. */
  typedef BCDThinPlateSplineKernelTransform Self;
  typedef itk::ThinPlateSplineKernelTransform<TScalarType,
                                              NDimension> Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** New macro for creation of through a Smart Pointer */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BCDThinPlateSplineKernelTransform,
               ThinPlateSplineKernelTransform);
  typename Superclass::Superclass::AMatrixType GetAMatrix()
  {
    return this->Superclass::Superclass::m_AMatrix;
  }

  typename Superclass::Superclass::BMatrixType GetBVector()
  {
    return this->Superclass::Superclass::m_BVector;
  }
};

#include "BRAINSLmkTransformCLP.h"

#include <fstream>
#include <vector>
#include <iostream>

/*
 * Description:
 *
 * This utility program estimates the affine transform to align the fixed landmarks
 * to the moving landmarks, and then generate the resampled moving image to the same
 * physical space as that of the reference image
 */

const   unsigned int ImageDimension = 3;
typedef short                                 PixelType;
typedef itk::Image<PixelType, ImageDimension> ImageType;
typedef std::vector<ImageType::PointType>     LandmarksVectorType;

LandmarksVectorType LoadLandmarks( std::string filename );

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  if( ( inputMovingLandmarks.compare( "" ) == 0 )
      && ( inputFixedLandmarks.compare( "" ) == 0 )
      && ( inputMovingVolume.compare( "" ) == 0 )
      && ( inputReferenceVolume.compare( "" ) == 0 ) )
    {
    itkGenericExceptionMacro(<< "Please set inputMovingLandmarks, inputFixedLandmarks, "
                             << "inputMovingVolume, and inputReferenceVolume.");
    }

  // typedefs
  typedef double CoordinateRepType;
  typedef itk::Point<
      CoordinateRepType, ImageDimension>           PointType;
  typedef std::vector<PointType>          PointArrayType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  typedef BCDThinPlateSplineKernelTransform<
      CoordinateRepType, ImageDimension>           TPSTransformType;
  typedef itk::AffineTransform<
      CoordinateRepType, ImageDimension>           AffineTransformType;
  typedef TPSTransformType::PointSetType PointSetType;
  typedef itk::TransformFileWriter       TransformWriterType;
  typedef PointSetType::Pointer          PointSetPointer;
  typedef PointSetType::PointIdentifier  PointIdType;
  typedef itk::ResampleImageFilter<
      ImageType, ImageType>                       ResamplerType;
  typedef itk::LinearInterpolateImageFunction<
      ImageType, double>                           InterpolatorType;

  // Read in landmarks
  PointSetType::Pointer sourceLandmarks = PointSetType::New();
  PointSetType::Pointer targetLandmarks = PointSetType::New();
    {
    PointSetType::PointsContainer::Pointer sourceLandmarkContainer =
      sourceLandmarks->GetPoints();
    PointSetType::PointsContainer::Pointer targetLandmarkContainer =
      targetLandmarks->GetPoints();
    PointIdType         id = itk::NumericTraits<PointIdType>::Zero;
    LandmarksVectorType targetLandmarksVec = LoadLandmarks( inputMovingLandmarks );
    LandmarksVectorType sourceLandmarksVec = LoadLandmarks( inputFixedLandmarks );

    // Sanity check
    if( targetLandmarksVec.size() != sourceLandmarksVec.size() )
      {
      std::cerr << "Different number of fixed and moving landmarks!" << std::endl;
      return EXIT_FAILURE;
      }

    unsigned int idx = 0;
    for( idx = 0; idx < sourceLandmarksVec.size(); ++idx )
      {
      sourceLandmarkContainer->InsertElement( id, sourceLandmarksVec[idx] );
      targetLandmarkContainer->InsertElement( id++, targetLandmarksVec[idx] );
      }
    }

  // Estimate affine transform
  AffineTransformType::Pointer affine = AffineTransformType::New();
    {
    TPSTransformType::Pointer tps = TPSTransformType::New();
    tps->SetSourceLandmarks( sourceLandmarks );
    tps->SetTargetLandmarks( targetLandmarks );
    tps->ComputeWMatrix();
    itk::Matrix<double, ImageDimension, ImageDimension> aMatrix( tps->GetAMatrix() );
    itk::Vector<double, ImageDimension>                 bVector;
    bVector.SetVnlVector( tps->GetBVector() );
    itk::Matrix<double, ImageDimension, ImageDimension> identity;
    identity.SetIdentity();
    affine->SetMatrix( aMatrix + identity );
    affine->SetOffset( bVector );
    }

  // Write output aligning transform
  if( outputAffineTransform.compare( "" ) != 0 )
    {
    TransformWriterType::Pointer writer = TransformWriterType::New();
    writer->SetInput( affine );
    writer->SetFileName( outputAffineTransform );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Cannot write the outputTransform file!" << std::endl;
      std::cerr << excep << std::endl;
      }
    }

  // Read in images
  ImageType::Pointer movingImage;
    {
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName( inputMovingVolume );
    try
      {
      reader->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }
    movingImage = reader->GetOutput();
    }

  ImageType::Pointer referenceImage;
    {
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName( inputReferenceVolume );
    try
      {
      reader->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }
    referenceImage = reader->GetOutput();
    }

  // Resample moving image
  ResamplerType::Pointer resampler = ResamplerType::New();
    {
    InterpolatorType::Pointer interpolator  = InterpolatorType::New();
    resampler->SetUseReferenceImage( true );
    resampler->SetInput( movingImage );
    resampler->SetReferenceImage( referenceImage );
    resampler->SetInterpolator( interpolator );
    resampler->SetTransform( affine );
    }

  // Write aligned image
  if( outputResampledVolume.compare( "" ) != 0 )
    {
    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput( resampler->GetOutput() );
    writer->SetFileName( outputResampledVolume );
    writer->SetUseCompression( true );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }
    }

  return EXIT_SUCCESS;
}

LandmarksVectorType
LoadLandmarks( std::string filename )
{
  LandmarksVectorType landmarks;
  std::string         line;
  std::ifstream       myfile( filename.c_str() );

  if( !myfile.is_open() )
    {
    itkGenericExceptionMacro(<< "Fatal error: Failed to load landmarks file. Program abort!");
    }
  while( getline(myfile, line) )
    {
    if( line.compare(0, 1, "#") != 0 )
      {
      unsigned int i;
      int          pos1 = line.find(',', 0);
      int          pos2;
      std::string  name = line.substr(0, pos1);
      if( name.compare("CM") == 0 )  // exclude CM
        {
        continue;
        }
      ImageType::PointType labelPos;
      for( i = 0; i < 3; ++i )
        {
        pos2 = line.find(',', pos1 + 1);
        labelPos[i] = atof( line.substr( pos1 + 1, pos2 - pos1 - 1 ).c_str() );
        if( i < 2 )
          {
          labelPos[i] *= -1; // RAS -> LPS
          }
        pos1 = pos2;
        }
      landmarks.push_back( labelPos );
      }
    }

  myfile.close();
  return landmarks;
}
