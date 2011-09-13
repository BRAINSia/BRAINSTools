#include "Utilities.h"

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

#include "itkBRAINSROIAutoImageFilter.h"

#include "GenericTransformImage.h"

template <typename TTransformType>
ProbabilityMapImageType::Pointer ItkTransformWarpImage(TTransformType transform,
                                                       ProbabilityMapImageType::Pointer & inputImage)
{
  typedef itk::ResampleImageFilter<ProbabilityMapImageType, ProbabilityMapImageType, double> ResampleFilterType;
  ResampleFilterType::Pointer         resampleFilter = ResampleFilterType::New();
  ProbabilityMapImageType::RegionType imageRegion = inputImage->GetBufferedRegion();
  resampleFilter->SetOutputParametersFromImage(inputImage);
  resampleFilter->SetInput(inputImage);
  resampleFilter->SetTransform(transform);

#if 0  // WindowedSinc is very slow, and almost certainly not necessary here.
  typedef itk::ConstantBoundaryCondition<ProbabilityMapImageType> BoundaryConditionType;
  const unsigned int WindowRadius = 5;
  typedef itk::Function::HammingWindowFunction<WindowRadius, double, double> WindowFunctionType;
  typedef itk::WindowedSincInterpolateImageFunction<
      ProbabilityMapImageType,
      WindowRadius,
      WindowFunctionType,
      BoundaryConditionType,
      double>    InterpolatorType;
#else
  typedef itk::LinearInterpolateImageFunction<ProbabilityMapImageType, double> InterpolatorType;
#endif

  InterpolatorType::Pointer interpolator  = InterpolatorType::New();
  resampleFilter->SetInterpolator(interpolator);
  std::cout << "Init resampling" << std::endl;
  try
    {
    resampleFilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception in Resampling." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }
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
    std::cerr << "Exception in CreateMITransformFile." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }
  return 0;
}

int CreateTransformFile(const std::string & MovingImageFilename,
                        const std::string & FixedImageFilename,
                        const std::string & OutputRegName,
                        const std::string & FixedBinaryImageFilename,  // ==
                                                                       // "NA"
                                                                       // if not
                                                                       // used
                        const std::string & MovingBinaryImageFilename, // ==
                                                                       // "NA"
                                                                       // if not
                                                                       // used
                                                                       // const
                                                                       // std::string
                                                                       // &
                                                                       // Command,
                                                                       // const
                                                                       // std::string
                                                                       // &
                                                                       // MovingAtlasFilename,
                                                                       // const
                                                                       // std::string
                                                                       // &
                                                                       // FixedAtlasFilename,
                                                                       // const
                                                                       // std::string
                                                                       // &
                                                                       // MovingAtlasBinaryFilename,
                                                                       // const
                                                                       // std::string
                                                                       // &
                                                                       // FixedAtlasBinaryFilename,
                                                                       //   bool
                                                                       // background,
                                                                       // bool
                                                                       // printRegistrationCommands,
                        bool verbose)
{
  std::cout << "* CreateTransformFile" << std::endl;
  std::cout << "MovingBinaryImageFilename :: "
            << MovingBinaryImageFilename << std::endl
            << "FixedBinaryImageFilename :: "
            << FixedBinaryImageFilename  << std::endl;

//
// START HERE MY TESTING WITH BRAINS HELPER ------------------------------------
// //
#if 1
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

  typedef itk::ImageFileReader<ProbabilityMapImageType> FixedVolumeReaderType;
  FixedVolumeReaderType::Pointer fixedVolumeReader = FixedVolumeReaderType::New();
  fixedVolumeReader->SetFileName( FixedImageFilename );
  try
    {
    fixedVolumeReader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception in fixedVolumeReader->Update();." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }

  ProbabilityMapImageType::Pointer fixedVolume = fixedVolumeReader->GetOutput();

  BSplineRegistrationHelper->SetFixedVolume( fixedVolume );

  // Set Moving Volume

  typedef itk::ImageFileReader<ProbabilityMapImageType> MovingVolumeReaderType;
  MovingVolumeReaderType::Pointer movingVolumeReader = MovingVolumeReaderType::New();
  movingVolumeReader->SetFileName( MovingImageFilename );
  try
    {
    movingVolumeReader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception in movingVolumeReader->Update();." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }

  ProbabilityMapImageType::Pointer movingVolume = movingVolumeReader->GetOutput();

  BSplineRegistrationHelper->SetMovingVolume( movingVolume );

  // Set Binary Images for BRAINSFit
  typedef itk::Image<unsigned char, 3> BinaryImageType;

  // - Fixed Image Binary Mask

  if( FixedBinaryImageFilename == "NA" ||
      FixedBinaryImageFilename == "na" ||
      FixedBinaryImageFilename == "" )
    {
    typedef itk::BRAINSROIAutoImageFilter<ProbabilityMapImageType,
                                          BinaryImageType> ROIAutoFilterType;

    ROIAutoFilterType::Pointer fixedVolumeROIFilter = ROIAutoFilterType::New();

    fixedVolumeROIFilter->SetInput( fixedVolume );
    fixedVolumeROIFilter->SetDilateSize(3);;
    try
      {
      fixedVolumeROIFilter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception in fixedVolumeROIFilter->Update();." << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      exit(-1);
      }

    BSplineRegistrationHelper->SetFixedBinaryVolume(
      fixedVolumeROIFilter->GetSpatialObjectROI() );
    // get brain region size
    }
  else
    {
    typedef itk::ImageFileReader<BinaryImageType> BinaryImageReaderType;
    BinaryImageReaderType::Pointer binaryFixedImageReader = BinaryImageReaderType::New();
    binaryFixedImageReader->SetFileName( FixedBinaryImageFilename );
    try
      {
      binaryFixedImageReader->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception in binaryFixedImageReader->Update();." << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      exit(-1);
      }
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
    typedef itk::BRAINSROIAutoImageFilter<ProbabilityMapImageType,
                                          BinaryImageType> ROIAutoFilterType;

    ROIAutoFilterType::Pointer movingVolumeROIFilter = ROIAutoFilterType::New();

    movingVolumeROIFilter->SetInput( movingVolume );
    movingVolumeROIFilter->SetDilateSize(1);;
    try
      {
      movingVolumeROIFilter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception in movingVolumeROIFilter->Update();." << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      exit(-1);
      }

    BSplineRegistrationHelper->SetMovingBinaryVolume(
      movingVolumeROIFilter->GetSpatialObjectROI() );
    }
  else
    {
    typedef itk::ImageFileReader<BinaryImageType> BinaryImageReaderType;
    BinaryImageReaderType::Pointer binaryMovingImageReader =
      BinaryImageReaderType::New();
    binaryMovingImageReader->SetFileName( MovingBinaryImageFilename );
    try
      {
      binaryMovingImageReader->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception in binaryMovingImageReader->Update();." << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      exit(-1);
      }

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

  try
    {
    BSplineRegistrationHelper->StartRegistration();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception in BSplineRegistrationHelper->StartRegistration();." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }
  catch( ... )
    {
    std::cerr << " Undefined Exception while BSplineRegistrationHelper->StartRegistration(); "
              << std::endl;
    exit(-1);
    }
  if( verbose > 0 )
    {
    std::cout << " - Write deformation " << std::endl
              << " :: " << OutputRegName
              << std::endl;
    }
  WriteTransformToDisk( BSplineRegistrationHelper->GetCurrentGenericTransform(),
                        OutputRegName );
  // Write out Transformed Output As Well
  // - EX. from GenericTransformImage.txx

  ProbabilityMapImageType::Pointer DeformedMovingImage;
  DeformedMovingImage = TransformResample<ProbabilityMapImageType,
                                          ProbabilityMapImageType>(
      movingVolume,
      fixedVolume,
      0.0F,
      GetInterpolatorFromString<ProbabilityMapImageType>("Linear"),
      BSplineRegistrationHelper->GetCurrentGenericTransform() );

  typedef itk::ImageFileWriter<ProbabilityMapImageType> DeformedVolumeWriterType;

  DeformedVolumeWriterType::Pointer deformedVolumeWriter = DeformedVolumeWriterType::New();

  deformedVolumeWriter->SetFileName( OutputRegName + "_output.nii.gz" );
  deformedVolumeWriter->SetInput( DeformedMovingImage );
  try
    {
    deformedVolumeWriter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception in deformedVolumeWriter->Update();." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }

#endif
// -------------------------------------- END HERE MY TESTING WITH BRAINS HELPER
// //

// START HERE  ----------------------------------------------------------------
// //
// Following part uses registration SCRIPT to generate transformation
//
//
#if 0

  // Mutualinformaiton code goes here
  //  std::cout << "Computing Reverse MI Registration." << std::endl;
  try
    {
    std::string command = Command;
    command += " ";
    command += MovingImageFilename;
    command += " ";
    command += FixedImageFilename;
    command += " ";
    command += OutputRegName;
    command += " ";
    command += MovingAtlasBinaryFilename;
    command += " ";
    command += FixedAtlasBinaryFilename;
    if( MovingAtlasFilename != "" )
      {
      if( FixedAtlasFilename == "" )
        {
        throw itk::ExceptionObject(__FILE__,
                                   __LINE__,
                                   "Need to specify both Moving and Fixed Landmarks",
                                   "CreateTransformFile");
        }
      command += " ";
      command += MovingAtlasFilename;
      command += " ";
      command += FixedAtlasFilename;
      }
    std::cout << "Running: " << command << std::endl;
    //    int status=system(command.c_str());
    itksysProcess *theProcess = itksysProcess_New();
    // Detach is for starting daemons, which we don't want to do...
    //    itksysProcess_SetOption(theProcess,itksysProcess_Option_Detach,1);
    itksysProcess_SetOption(theProcess, itksysProcess_Option_HideWindow, 1);
    itksysProcess_SetOption(theProcess, itksysProcess_Option_Verbatim, 1);
    itksysProcess_SetPipeFile(theProcess, itksysProcess_Pipe_STDOUT, "/dev/null");
    itksysProcess_SetPipeFile(theProcess, itksysProcess_Pipe_STDERR, "/dev/null");
    //    itksysProcess_SetPipeShared(theProcess,itksysProcess_Pipe_STDOUT,1);
    //    itksysProcess_SetPipeShared(theProcess,itksysProcess_Pipe_STDERR,1);
    const char *argv[2];
    argv[0] = command.c_str();
    argv[1] = 0;
    itksysProcess_SetCommand(theProcess, argv);
    itksysProcess_Execute(theProcess);
    *process_handle = theProcess;
    return 0;
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception in Transformation Creation." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }
#endif
// ---------------------------------------------------------------------- END
// HERE//
  return 1;
}

#if 0
// THIS IS NO LONGER NEEDED!!
bool CHECK_CORONAL(ProbabilityMapImageType::DirectionType Dir)
{
  if( Dir !=
      itk::SpatialOrientationAdapter().ToDirectionCosines(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP) )
    {
    std::cout << "Image Directions are not in RIP orientation "
              << std::endl
              << "\nRequire\n"
              <<  itk::SpatialOrientationAdapter().ToDirectionCosines(
      itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP)
              << std::endl;
    return false;
    }
  return true;
}

#endif

void CreateNewImageFromTemplate(ProbabilityMapImageType::Pointer & PointerToOutputImage,
                                const ProbabilityMapImageType::Pointer & PreInitializedImage)
{
  ProbabilityMapImageType::RegionType region;

  PointerToOutputImage = ProbabilityMapImageType::New();
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
  itk::ImageRegionIterator<ProbabilityMapImageType> bbri( PointerToOutputImage,
                                                          PointerToOutputImage->GetLargestPossibleRegion() );
  bbri = bbri.Begin();
  while( !bbri.IsAtEnd() )
    {
    // Zeroing voxel signal intensity values
    bbri.Set(itk::NumericTraits<ProbabilityMapImageType::PixelType>::Zero);
    ++bbri;
    }
}

void CreateNewFloatImageFromTemplate(RealImageType::Pointer & PointerToOutputImage,
                                     const ProbabilityMapImageType::Pointer & PreInitializedImage)
{
  RealImageType::RegionType region;

  PointerToOutputImage = RealImageType::New();
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
  itk::ImageRegionIterator<RealImageType> bbri( PointerToOutputImage,
                                                PointerToOutputImage->GetLargestPossibleRegion() );
  bbri = bbri.Begin();
  while( !bbri.IsAtEnd() )
    {
    // Zeroing voxel signal intensity values
    bbri.Set(itk::NumericTraits<RealImageType::PixelType>::Zero);
    ++bbri;
    }
}

// REGINA:: removed #if 0

void ReadDeformationField(TDeformationField::Pointer & DeformationField,
                          const std::string & DeformationFilename)
{
  std::cout << "Reading Deformation Field " << DeformationFilename.c_str() << std::endl;

  try
    {
    itk::ImageFileReader<TDeformationField>::Pointer DeformationReader =
      itk::ImageFileReader<TDeformationField>::New();
    DeformationReader->SetFileName( DeformationFilename.c_str() );
    DeformationReader->Update();
    DeformationField = DeformationReader->GetOutput();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception in reading deformation field " << DeformationFilename.c_str() << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }
  catch( std::exception & error )
    {
    std::cout << "Memory Error near " << __FILE__ <<  " " << __LINE__ << ": " << error.what() << std::endl;
    exit(-1);
    }
  catch( ... )
    {
    std::cout << "Unknown error near "  << __FILE__ <<  " " << __LINE__ << std::endl;
    exit(-1);
    }
  // REGINA:: removed #if 0 //HACK:  For right now, everything must be in RIP
  // orientaiton to work.
  std::cout << "Deformation Field Orientation \n" << DeformationField->GetDirection() << std::endl << std::endl;
}

/*
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
*/

void FillGradProfile(std::vector<float>::iterator & fi,
                     const std::map<std::string, ProbabilityMapImageType::Pointer> MapOfImages,
                     std::map<std::string, ImageLinearInterpolatorType::Pointer>  MapOfImageInterpolators,
                     //                     const
                     // ProbabilityMapImageType::Pointer DeformedProbMap,
                     const itk::Image<itk::CovariantVector<float, 3>, 3>::Pointer ProbMapGradient,
                     const ProbabilityMapImageType::IndexType & CurrentIndex,
                     const int ProfileExtent/*,
                                              * const int IrisSize   --> We are
                                              * not using iris size anymore!*/)
{
  for( std::map<std::string, ProbabilityMapImageType::Pointer>::const_iterator imapi = MapOfImages.begin();
       imapi != MapOfImages.end();
       imapi++ )
    {
    const ProbabilityMapImageType::SpacingType ImageSpacing = imapi->second->GetSpacing();
    const float                                deltax = ProbMapGradient->GetPixel(CurrentIndex)[0];
    const float                                deltay = ProbMapGradient->GetPixel(CurrentIndex)[1];
    const float                                deltaz = ProbMapGradient->GetPixel(CurrentIndex)[2];

    const float GradientVectorLength = vcl_sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
    const float InvGradientVectorLength = ( GradientVectorLength > 0.0F ) ? 1.0 / GradientVectorLength : 1;
    const float unitdeltax = deltax * InvGradientVectorLength;
    const float unitdeltay = deltay * InvGradientVectorLength;
    const float unitdeltaz = deltaz * InvGradientVectorLength;
      {
      itk::Point<float, 3> BasePhysicalPoint;
      imapi->second->TransformIndexToPhysicalPoint(CurrentIndex, BasePhysicalPoint);
      for( float i = -ProfileExtent; i <= ProfileExtent; i += 1.0F )
        {
        itk::Point<float, 3> CurrProfilePoint = BasePhysicalPoint;
        CurrProfilePoint[0] = BasePhysicalPoint[0] + i * unitdeltax;
        CurrProfilePoint[1] = BasePhysicalPoint[1] + i * unitdeltay;
        CurrProfilePoint[2] = BasePhysicalPoint[2] + i * unitdeltaz;
        itk::ContinuousIndex<float, 3> ContinuousIndexProfilePoint;
        const bool                     inBuffer =
          imapi->second
          ->TransformPhysicalPointToContinuousIndex(CurrProfilePoint, ContinuousIndexProfilePoint);
        // From the imageinterpolators, Get the continuous index value.
        if( inBuffer )  //
        /*
          * MapOfImageInterpolators[imapi->first]->IsInsideBuffer(ContinuousIndexProfilePoint))
          */
          {
          *fi =  static_cast<float>( MapOfImageInterpolators[imapi->first]
                                     ->EvaluateAtContinuousIndex(ContinuousIndexProfilePoint) );

          //        std::cout<<ContinuousIndexProfilePoint<<" :
          // ["<<*fi<<","<<MapOfImageInterpolators[imapi->first]
          //            ->EvaluateAtContinuousIndex(ContinuousIndexProfilePoint)
          // <<"],";
          }
        else
          {
          std::cerr << "THIS IS HIGHLY UNLIKELY TO HAPPEN, AND NOT ALLOWED."
                    << " The error is occuring at location " << CurrentIndex
                    << " while accessing point " << ContinuousIndexProfilePoint << " with unit delta "
                    << unitdeltax << " " << unitdeltay << " " << unitdeltaz << std::endl;
          std::cerr << "Trying to access point :"
                    << ContinuousIndexProfilePoint << " in image where physical space is from \n"
                    << imapi->second->GetOrigin() << " to [";
          for( int q = 0; q < 3; q++ )
            {
            std::cerr << imapi->second->GetOrigin()[q]
              + imapi->second->GetLargestPossibleRegion().GetSize()[q]
              * imapi->second->GetSpacing()[q] << ",";
            }
          std::cerr << "]" << std::endl;
          exit(-1);
          *fi = 0.0F;
          }

        fi++;
        }
      }
    }
  return;
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

void AddInputVector
  (std::vector<float> & inputvector,
  // ProbabilityMapImageType::Pointer ProbMapImage,
  itk::Image<itk::CovariantVector<float, 3>, 3>::Pointer ProbMapGradient,
  RealImageType::Pointer RhoMapImage,
  RealImageType::Pointer PhiMapImage,
  RealImageType::Pointer ThetaMapImage,
  // const DataSet::TypeVector &ImageTypeList,
  const std::map<std::string, ProbabilityMapImageType::Pointer> & MapOfImages,
  std::map<std::string, ImageLinearInterpolatorType::Pointer>  & MapOfImageInterpolators,
  const ProbabilityMapImageType::Pointer DeformedProbabilityMap[],
  const int NumberOfProbabilityMaps,
  const ProbabilityMapImageType::IndexType & CurrentIndex,
  const int GradientProfileSize)
{
  // std::cout<<"Add input Vector...."<<std::endl;
  std::vector<float>::iterator fi = inputvector.begin();
  // std::cout<<"EMbrased part start here.....:"<<std::endl;
    {
    *fi = 1.0F; // Hardcode to 1.0 just to eliminate it's effect
    // static_cast<float>(ProbMapImage->GetPixel(CurrentIndex));
    // *fi=static_cast<float>(ProbMapImage->GetPixel(CurrentIndex));
    // std::cout<<"\n"<<*fi<<",";
    // fi++;
    /*========================================================================//
      * Add Sphericcal Coord System
      * =========================================================================*/
    *fi = static_cast<float>( RhoMapImage->GetPixel(CurrentIndex) );
    // std::cout<<std::endl<<*fi<<",";
    ++fi;
    *fi = static_cast<float>( PhiMapImage->GetPixel(CurrentIndex) );
    // std::cout<<*fi<<",";
    ++fi;
    *fi = static_cast<float>( ThetaMapImage->GetPixel(CurrentIndex) );
    // std::cout<<*fi<<",";
    ++fi;
    /*========================================================================//
      * Add Probability Maps....
      * =========================================================================*/
    for( int i = 0; i < NumberOfProbabilityMaps; i++ )
      {
      float probabilityValue = DeformedProbabilityMap[i]->GetPixel(CurrentIndex);

      if( probabilityValue > 0.0 )
        {
        *fi = 1.0F;
        // std::cout <<"Prob Value:::: " << probabilityValue << "== " << *fi <<"
        // ";
        }
      else
        {
        *fi = 0.0F;
        }
      ++fi;
      }
    // std::cout<<std::endl;
    }

  /*========================================================================//
    * Find Gradient Value....
    * =========================================================================*/
  FillGradProfile(fi,                       // Index + rho phi theta
                  MapOfImages,
                  MapOfImageInterpolators,
                  // ProbMapImage,
                  ProbMapGradient,
                  CurrentIndex,
                  GradientProfileSize/*,
                                       * IrisSize--> We are not using iris size
                                       * any more.  */);
}

void DefineBoundingBox(const ProbabilityMapImageType::Pointer image,
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
    ProbabilityMapImageType::RegionType EntireRegion;
    EntireRegion.SetSize( image->GetLargestPossibleRegion().GetSize() );
    EntireRegion.SetIndex( image->GetLargestPossibleRegion().GetIndex() );
    itk::ImageRegionIterator<ProbabilityMapImageType> ri(image, EntireRegion);
    while( !ri.IsAtEnd() )
      {
      if( ri.Get() > 0 )
        {
        const ProbabilityMapImageType::IndexType index = ri.GetIndex();
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

TDeformationField::Pointer
RegisterLandmarksToDeformationField(const std::string & InputImageFilename,
                                    const std::string & InputAtlasFilename,
                                    const std::string & TemplateAtlasFilename,
                                    const std::string & TemplateImageFilename)
{
  return itkUtil::RegisterLandmarksToDeformationField
         <TDeformationField, ProbabilityMapImageType>(InputImageFilename,
                                                      InputAtlasFilename,
                                                      TemplateAtlasFilename,
                                                      TemplateImageFilename);
}

#endif
