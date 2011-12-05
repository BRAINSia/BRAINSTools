#include "Utilities.h"
#include <NeuralParams.h>
#include "ANNParams.h"
#include "itksys/SystemTools.hxx"
#include "itkSpatialOrientationAdapter.h"
#include "itkResampleImageFilter.h"
#include "itkBSplineDeformableTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkQuaternionRigidTransform.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"
#include "itkAffineTransform.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "itkEuler3DTransform.h"
#include "itkImageMaskSpatialObject.h"
#include "itkMedianImageFilter.h"

#include "itkBRAINSROIAutoImageFilter.h"

#include "GenericTransformImage.h"
#include <sstream>

template <typename TTransformType>
InternalImageType::Pointer ItkTransformWarpImage(TTransformType transform,
                                                 InternalImageType::Pointer & inputImage)
{
  typedef itk::ResampleImageFilter<InternalImageType, InternalImageType, double> ResampleFilterType;
  ResampleFilterType::Pointer   resampleFilter = ResampleFilterType::New();
  InternalImageType::RegionType imageRegion = inputImage->GetBufferedRegion();
  resampleFilter->SetOutputParametersFromImage(inputImage);
  resampleFilter->SetInput(inputImage);
  resampleFilter->SetTransform(transform);

  typedef itk::LinearInterpolateImageFunction<InternalImageType, double> InterpolatorType;

  InterpolatorType::Pointer interpolator  = InterpolatorType::New();
  resampleFilter->SetInterpolator(interpolator);
  std::cout << "Init resampling" << std::endl;
  resampleFilter->Update();
  std::cout << "Resampled Image" << std::endl;
  return resampleFilter->GetOutput();
}

// REGINA:: removed #if 0

int CreateMITransformFile(const std::string & MovingImageFilename,
                          const std::string & FixedImageFilename,
                          const std::string & MovingImageToFixedImageRegistrationFilename,
                          const std::string & ParamName)
{
  // Mutualinformaiton code goes here
  std::cout << "Computing Reverse MI Registration." << std::endl;

  try
    {
    std::string command = ParamName + " " + MovingImageFilename + " " + FixedImageFilename + " "
      + MovingImageToFixedImageRegistrationFilename;
    std::cout << "Running: " << command << std::endl;
    int status = system( command.c_str() );
    return status;
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << e << std::endl;
    throw e;
    }
  return 0;
}

int CreateTransformFile(const std::string & MovingImageFilename,
                        const std::string & FixedImageFilename,
                        const std::string & OutputRegName,
                        const std::string & FixedBinaryImageFilename,
                        const std::string & MovingBinaryImageFilename,
                        const int   roiAutoDilateSize,
                        bool verbose)
{
  std::cout << "* CreateTransformFile" << std::endl;
  std::cout << "MovingBinaryImageFilename :: "
            << MovingBinaryImageFilename << std::endl
            << "FixedBinaryImageFilename :: "
            << FixedBinaryImageFilename  << std::endl;

  // Create Helper Class of BRAINSFit for BSpline Registraion.
  // CAUTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Check if input images are all 'float'
  // BRAINSFit takes only float type images!
  // For now, BRAINSCut only uses float type, so it is safe!!!!!!

  typedef itk::BRAINSFitHelper BSplineRegistrationHelperType;
  BSplineRegistrationHelperType::Pointer BSplineRegistrationHelper =
    BSplineRegistrationHelperType::New();

  // Set BSpline Trainsformation Type
  std::vector<std::string> transformType( 4 );

  // transformType[0]="Rigid";
  transformType[0] = "ScaleVersor3D";
  transformType[1] = "ScaleSkewVersor3D";
  transformType[2] = "Affine";
  transformType[3] = "BSpline";

  BSplineRegistrationHelper->SetTransformType( transformType );

  // Set Displacement Field Size
  BSplineRegistrationHelper->SetMaxBSplineDisplacement(4);

  // BSpline Grid Size
  std::vector<int> splineGridSize(3);

  splineGridSize[0] = 28;
  splineGridSize[1] = 20;
  splineGridSize[2] = 24;

  BSplineRegistrationHelper->SetSplineGridSize( splineGridSize );

  // Number of Iterations
  std::vector<int> numberOfIterations(4);

  numberOfIterations[0] = 1500;
  numberOfIterations[1] = 1500;
  numberOfIterations[2] = 1500;
  numberOfIterations[3] = 1500;
  // numberOfIterations[4]=1500;

  BSplineRegistrationHelper->SetNumberOfIterations( numberOfIterations );

  // Set Minimum Step Length
  std::vector<double> minimumStepLength(4);

  minimumStepLength[0] = 0.005;
  minimumStepLength[1] = 0.005;
  minimumStepLength[2] = 0.005;
  minimumStepLength[3] = 0.005;
  // minimumStepLength[4]=0.005;

  BSplineRegistrationHelper->SetMinimumStepLength( minimumStepLength );

  // Set Histogram Matching Options
  // BSplineRegistrationHelper->SetNumberOfHistogramBins(250);
  // TODO: decide what to do here too much freedom?
  // BSplineRegistrationHelper->SetNumberOfMatchPoints(10);
  // TODO:: HISTOGRAMMATCHING OPTION is now Disabled!
  BSplineRegistrationHelper->SetHistogramMatch(false);

  // Set Fixed Volume

  typedef itk::ImageFileReader<InternalImageType> FixedVolumeReaderType;
  FixedVolumeReaderType::Pointer fixedVolumeReader = FixedVolumeReaderType::New();
  fixedVolumeReader->SetFileName( FixedImageFilename );

  fixedVolumeReader->Update();

  InternalImageType::Pointer fixedVolume = fixedVolumeReader->GetOutput();

  BSplineRegistrationHelper->SetFixedVolume( fixedVolume );

  // Set Moving Volume

  typedef itk::ImageFileReader<InternalImageType> MovingVolumeReaderType;
  MovingVolumeReaderType::Pointer movingVolumeReader = MovingVolumeReaderType::New();
  movingVolumeReader->SetFileName( MovingImageFilename );

  movingVolumeReader->Update();

  InternalImageType::Pointer movingVolume = movingVolumeReader->GetOutput();

  BSplineRegistrationHelper->SetMovingVolume( movingVolume );

  typedef itk::Image<unsigned char, 3> BinaryImageType;

  // - Fixed Image Binary Mask

  if( FixedBinaryImageFilename == "NA" ||
      FixedBinaryImageFilename == "na" ||
      FixedBinaryImageFilename == "" )
    {
    typedef itk::BRAINSROIAutoImageFilter<InternalImageType,
                                          BinaryImageType> ROIAutoFilterType;

    ROIAutoFilterType::Pointer fixedVolumeROIFilter = ROIAutoFilterType::New();

    fixedVolumeROIFilter->SetInput( fixedVolume );
    fixedVolumeROIFilter->SetDilateSize(roiAutoDilateSize);
    fixedVolumeROIFilter->Update();

    BSplineRegistrationHelper->SetFixedBinaryVolume(
      fixedVolumeROIFilter->GetSpatialObjectROI() );
    // get brain region size
    }
  else
    {
    typedef itk::ImageFileReader<BinaryImageType> BinaryImageReaderType;
    BinaryImageReaderType::Pointer binaryFixedImageReader = BinaryImageReaderType::New();
    binaryFixedImageReader->SetFileName( FixedBinaryImageFilename );
    binaryFixedImageReader->Update();

    typedef itk::ImageMaskSpatialObject<3> binarySpatialObjectType;
    binarySpatialObjectType::Pointer binaryFixedObject
      = binarySpatialObjectType::New();
    binaryFixedObject->SetImage( binaryFixedImageReader->GetOutput() );

    BSplineRegistrationHelper->SetFixedBinaryVolume( binaryFixedObject );
    }

  // - Moving Image Binary Mask

  if( MovingBinaryImageFilename == "NA" ||
      MovingBinaryImageFilename == "na" ||
      MovingBinaryImageFilename == "" )
    {
    typedef itk::BRAINSROIAutoImageFilter<InternalImageType,
                                          BinaryImageType> ROIAutoFilterType;

    ROIAutoFilterType::Pointer movingVolumeROIFilter = ROIAutoFilterType::New();

    movingVolumeROIFilter->SetInput( movingVolume );
    movingVolumeROIFilter->SetDilateSize(roiAutoDilateSize);
    movingVolumeROIFilter->Update();

    BSplineRegistrationHelper->SetMovingBinaryVolume(
      movingVolumeROIFilter->GetSpatialObjectROI() );
    }
  else
    {
    typedef itk::ImageFileReader<BinaryImageType> BinaryImageReaderType;
    BinaryImageReaderType::Pointer binaryMovingImageReader =
      BinaryImageReaderType::New();
    binaryMovingImageReader->SetFileName( MovingBinaryImageFilename );
    binaryMovingImageReader->Update();

    std::cout << " Set Binary Image " << std::endl;
    typedef itk::ImageMaskSpatialObject<3> binarySpatialObjectType;
    binarySpatialObjectType::Pointer binaryMovingObject
      = binarySpatialObjectType::New();
    binaryMovingObject->SetImage( binaryMovingImageReader->GetOutput() );

    BSplineRegistrationHelper->SetMovingBinaryVolume( binaryMovingObject );
    }

  // Set Other Options

  // TODO :: Following has to be user input
  const unsigned int numberOfSamples = 100000;
  const unsigned int maxBSplineDisplacement = 7;
  const double       maskInferiorCutOffFromCenter = 65.0;
  const double       translationScale = 1000.0;
  const double       reproportionalScale = 1.0;
  const double       skewScale = 1.0;

  BSplineRegistrationHelper->SetNumberOfSamples( numberOfSamples );
  BSplineRegistrationHelper->SetMaxBSplineDisplacement( maxBSplineDisplacement );
  BSplineRegistrationHelper->SetInitializeTransformMode( "useCenterOfHeadAlign" );
  BSplineRegistrationHelper->SetMaskInferiorCutOffFromCenter(
    maskInferiorCutOffFromCenter );
  // BSplineRegistrationHelper->SetOutputTransform( OutputRegName );
  BSplineRegistrationHelper->SetTranslationScale( translationScale );
  BSplineRegistrationHelper->SetReproportionScale( reproportionalScale );
  BSplineRegistrationHelper->SetSkewScale( skewScale );

  //  Start Registration

  // TODO: is this line really print before start registration????

  if( verbose > 0 )
    {
    BSplineRegistrationHelper->PrintCommandLine(true, "BSplineRegistrationHelper");
    }

  BSplineRegistrationHelper->StartRegistration();

  if( verbose > 0 )
    {
    std::cout << " - Write deformation " << std::endl
              << " :: " << OutputRegName
              << std::endl;
    }
  WriteTransformToDisk( BSplineRegistrationHelper->GetCurrentGenericTransform(),
                        OutputRegName );
  // Write out Transformed Output As Well
  // - EX. from GenericTransformImage.hxx

  InternalImageType::Pointer DeformedMovingImage;
  DeformedMovingImage = TransformResample<InternalImageType,
                                          InternalImageType>(
      movingVolume,
      fixedVolume,
      0.0F,
      GetInterpolatorFromString<InternalImageType>("Linear"),
      BSplineRegistrationHelper->GetCurrentGenericTransform() );

  typedef itk::ImageFileWriter<InternalImageType> DeformedVolumeWriterType;

  DeformedVolumeWriterType::Pointer deformedVolumeWriter = DeformedVolumeWriterType::New();

  deformedVolumeWriter->SetFileName( OutputRegName + "_output.nii.gz" );
  deformedVolumeWriter->SetInput( DeformedMovingImage );
  deformedVolumeWriter->Update();
// -------------------------------------- END HERE MY TESTING WITH BRAINS HELPER

  return 1;
}

void CreateNewImageFromTemplate(InternalImageType::Pointer & PointerToOutputImage,
                                const InternalImageType::Pointer & PreInitializedImage)
{
  InternalImageType::RegionType region;

  PointerToOutputImage = InternalImageType::New();
  region.SetSize( PreInitializedImage->GetLargestPossibleRegion().GetSize() );
  region.SetIndex( PreInitializedImage->GetLargestPossibleRegion().GetIndex() );
  PointerToOutputImage->SetLargestPossibleRegion(region);
  PointerToOutputImage->SetBufferedRegion(region);
  PointerToOutputImage->SetRequestedRegion(region);
  PointerToOutputImage->CopyInformation(PreInitializedImage);
  PointerToOutputImage->Allocate();
  PointerToOutputImage->FillBuffer(0.0);
  PointerToOutputImage->SetDirection( PreInitializedImage->GetDirection() );
  // NO LONGER NEEDED CHECK_CORONAL(PointerToOutputImage->GetDirection());
  /*
    * PointerToOutputImage->SetMetaDataDictionary(PreInitializedImage->GetMetaDataDictionary());
    */
  itk::ImageRegionIterator<InternalImageType> bbri( PointerToOutputImage,
                                                    PointerToOutputImage->GetLargestPossibleRegion() );
  bbri = bbri.Begin();
  while( !bbri.IsAtEnd() )
    {
    // Zeroing voxel signal intensity values
    bbri.Set(itk::NumericTraits<InternalImageType::PixelType>::Zero);
    ++bbri;
    }
}

void CreateNewFloatImageFromTemplate(InternalImageType::Pointer & PointerToOutputImage,
                                     const InternalImageType::Pointer & PreInitializedImage)
{
  InternalImageType::RegionType region;

  PointerToOutputImage = InternalImageType::New();
  region.SetSize( PreInitializedImage->GetLargestPossibleRegion().GetSize() );
  region.SetIndex( PreInitializedImage->GetLargestPossibleRegion().GetIndex() );
  PointerToOutputImage->SetLargestPossibleRegion(region);
  PointerToOutputImage->SetBufferedRegion(region);
  PointerToOutputImage->SetRequestedRegion(region);
  PointerToOutputImage->CopyInformation(PreInitializedImage);
  PointerToOutputImage->Allocate();
  PointerToOutputImage->FillBuffer(0.0);
  // NO LONGER NEEDED CHECK_CORONAL(PointerToOutputImage->GetDirection());
  PointerToOutputImage->SetDirection( PreInitializedImage->GetDirection() );
  PointerToOutputImage->SetMetaDataDictionary( PreInitializedImage->GetMetaDataDictionary() );
  itk::ImageRegionIterator<InternalImageType> bbri( PointerToOutputImage,
                                                    PointerToOutputImage->GetLargestPossibleRegion() );
  bbri = bbri.Begin();
  while( !bbri.IsAtEnd() )
    {
    // Zeroing voxel signal intensity values
    bbri.Set(itk::NumericTraits<InternalImageType::PixelType>::Zero);
    ++bbri;
    }
}

// REGINA:: removed #if 0

void ReadDisplacementField(TDisplacementField::Pointer & DisplacementField,
                           const std::string & DisplacementFilename)
{
  std::cout << "Reading Displacement Field " << DisplacementFilename.c_str() << std::endl;

  try
    {
    itk::ImageFileReader<TDisplacementField>::Pointer DeformationReader =
      itk::ImageFileReader<TDisplacementField>::New();
    DeformationReader->SetFileName( DisplacementFilename.c_str() );
    DeformationReader->Update();
    DisplacementField = DeformationReader->GetOutput();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << e << std::endl;
    throw e;
    }

  // REGINA:: removed #if 0 //HACK:  For right now, everything must be in RIP
  // orientaiton to work.
  std::cout << "Displacement Field Orientation \n" << DisplacementField->GetDirection() << std::endl << std::endl;
}

bool CheckIndex(const int SizeX,
                const int SizeY,
                const int SizeZ,
                const int X,
                const int Y,
                const int Z)
{
  if( ( X < 0 ) || ( X >= SizeX ) || ( Y < 0 ) || ( Y >= SizeY ) || ( Z < 0 ) || ( Z >= SizeZ ) )
    {
    return false;
    }
  return true;
}

bool CheckIndexFloat(const int SizeX,
                     const int SizeY,
                     const int SizeZ,
                     const float X,
                     const float Y,
                     const float Z)
{
  if( ( X < 0 ) || ( X >= SizeX ) || ( Y < 0 ) || ( Y >= SizeY ) || ( Z < 0 ) || ( Z >= SizeZ ) )
    {
    return false;
    }
  return true;
}

/*BUG:REGINA
  * PHI and THETA should be switched in their computing way to follow standard
  * definition of spherical coordinate system
  */
void XYZToSpherical(const itk::Point<float, 3> & LocationWithOriginAtCenterOfImage,
                    float & rho, float & phi, float & theta)
{
  /*Rho*/
#define _SQR(a) ( ( a ) * ( a ) )
  rho = static_cast<float>
    ( vcl_sqrt( _SQR(LocationWithOriginAtCenterOfImage[0])
                + _SQR(LocationWithOriginAtCenterOfImage[1])
                + _SQR(LocationWithOriginAtCenterOfImage[2]) ) );
#undef _SQR
  /*Phi*/
  phi = 0.0F;
  if( LocationWithOriginAtCenterOfImage[0] < 0 )
    {
    phi = vcl_atan2(-LocationWithOriginAtCenterOfImage[0], LocationWithOriginAtCenterOfImage[1]);
    }
  else
    {
    phi = vcl_atan2(LocationWithOriginAtCenterOfImage[0], LocationWithOriginAtCenterOfImage[1]);
    }
  /*Theta*/
  theta = 0.0F;
  if( LocationWithOriginAtCenterOfImage[2] < 0 )
    {
    theta = vcl_atan2(-LocationWithOriginAtCenterOfImage[2], LocationWithOriginAtCenterOfImage[1]);
    }
  else
    {
    theta = vcl_atan2(LocationWithOriginAtCenterOfImage[2], LocationWithOriginAtCenterOfImage[1]);
    }

  //  theta = vcl_acos(LocationWithOriginAtCenterOfImage[2]/rho);

  rho = rho / 128.0F;  // The largest brain ever will always fit in a sphere
                       // with radius of 128MM centered at the AC point
  phi = phi / (vnl_math::pi);
  theta = theta / (vnl_math::pi);
}

unsigned int InputVectorSizeRequirement(
  const unsigned int NumberOfProbabilityMaps,
  const unsigned int NumberOfImageTypes,
  const unsigned int GradientProfileSize
  )
{
  const unsigned int NumImageCoords = 3;
  const unsigned int ivsize = NumImageCoords + NumberOfProbabilityMaps
    + ( ( GradientProfileSize * 2 ) + 1 ) * NumberOfImageTypes;

  // const unsigned int ivsize=
  // NumImageCoords+((GradientProfileSize*2)+1)*NumberOfImageTypes;
  std::cout << "NumImageCoords: " << NumImageCoords << std::endl;
  std::cout << "NumberOfProbabilityMaps(Structure of Interest): " << NumberOfProbabilityMaps << std::endl;
  std::cout << "GradientProfileSize: " << GradientProfileSize << std::endl;
  std::cout << "NumberOfImageTypes: " << NumberOfImageTypes << std::endl;
  std::cout << "ivsize: " << ivsize << std::endl;
  std::cout << "ivsize= NumImageCoords+NumberOfProbabilityMaps+((GradientProfileSize*2)+1)*NumberOfImageTypes"
            << std::endl;

  return ivsize;
}

void DefineBoundingBox(const InternalImageType::Pointer image,
                       itk::Index<3> & min,
                       itk::Index<3> & max)
{
  const int NX = image->GetLargestPossibleRegion().GetSize()[0];
  const int NY = image->GetLargestPossibleRegion().GetSize()[1];
  const int NZ = image->GetLargestPossibleRegion().GetSize()[2];

  max[0] = 0;
  max[1] = 0;
  max[2] = 0;
  min[0] = NX;
  min[1] = NY;
  min[2] = NZ;
    {
    InternalImageType::RegionType EntireRegion;
    EntireRegion.SetSize( image->GetLargestPossibleRegion().GetSize() );
    EntireRegion.SetIndex( image->GetLargestPossibleRegion().GetIndex() );
    itk::ImageRegionIterator<InternalImageType> ri(image, EntireRegion);
    while( !ri.IsAtEnd() )
      {
      if( ri.Get() > 0 )
        {
        const InternalImageType::IndexType index = ri.GetIndex();
        min[0] = ( min[0] < index[0] ) ? min[0] : index[0];
        min[1] = ( min[1] < index[1] ) ? min[1] : index[1];
        min[2] = ( min[2] < index[2] ) ? min[2] : index[2];
        max[0] = ( max[0] > index[0] ) ? max[0] : index[0];
        max[1] = ( max[1] > index[1] ) ? max[1] : index[1];
        max[2] = ( max[2] > index[2] ) ? max[2] : index[2];
        }
      ++ri;
      }
    }
}

#if 0

TDisplacementField::Pointer
RegisterLandmarksToDeformationField(const std::string & InputImageFilename,
                                    const std::string & InputAtlasFilename,
                                    const std::string & TemplateAtlasFilename,
                                    const std::string & TemplateImageFilename)
{
  return itkUtil::RegisterLandmarksToDeformationField
         <TDisplacementField, InternalImageType>(InputImageFilename,
                                                 InputAtlasFilename,
                                                 TemplateAtlasFilename,
                                                 TemplateImageFilename);
}

#endif
