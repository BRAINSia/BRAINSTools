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
    throw;
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
  // - EX. from GenericTransformImage.txx

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
    throw;
    }
#endif
// ---------------------------------------------------------------------- END
// HERE//
  return 1;
}

#if 0
// THIS IS NO LONGER NEEDED!!
bool CHECK_CORONAL(InternalImageType::DirectionType Dir)
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
    throw;
    }

  // REGINA:: removed #if 0 //HACK:  For right now, everything must be in RIP
  // orientaiton to work.
  std::cout << "Deformation Field Orientation \n" << DeformationField->GetDirection() << std::endl << std::endl;
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

void FillGradProfile(std::vector<neural_scalar_type>::iterator & fi,
                     const std::map<std::string, InternalImageType::Pointer> & MapOfImages,
                     std::map<std::string, neural_scalar_type>& MapOfLocalMax,
                     std::map<std::string, neural_scalar_type>& MapOfLocalMin,
                     const itk::Image<itk::CovariantVector<float, 3>, 3>::Pointer ProbMapGradient,
                     const InternalImageType::IndexType & CurrentIndex,
                     const int ProfileExtent/*,
                                              * const int IrisSize   --> We are
                                              * not using iris size anymore!*/)
{
  // std::cout<<__LINE__<<" :: "<<__FILE__<< " :: PrfileExtent :: " <<ProfileExtent<<std::endl;
  ImageLinearInterpolatorType::Pointer Interpolator =
    ImageLinearInterpolatorType::New();

  // std::cout<<__LINE__<<" :: "<<__FILE__<<std::endl;
  for( std::map<std::string, InternalImageType::Pointer>::const_iterator imapi = MapOfImages.begin();
       imapi != MapOfImages.end();
       imapi++ )
    {
    // std::cout<<__LINE__<<" :: "<<__FILE__<<"::"<<imapi->second<<std::endl;
    Interpolator->SetInputImage( imapi->second );
    // std::cout<<__LINE__<<" :: "<<__FILE__<<std::endl;
    const InternalImageType::SpacingType ImageSpacing = imapi->second->GetSpacing();
    const float                          deltax = ProbMapGradient->GetPixel(CurrentIndex)[0];
    const float                          deltay = ProbMapGradient->GetPixel(CurrentIndex)[1];
    const float                          deltaz = ProbMapGradient->GetPixel(CurrentIndex)[2];

    const float GradientVectorLength = vcl_sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
    const float InvGradientVectorLength = ( GradientVectorLength > 0.0F ) ? 1.0 / GradientVectorLength : 1;
    const float unitdeltax = deltax * InvGradientVectorLength;
    const float unitdeltay = deltay * InvGradientVectorLength;
    const float unitdeltaz = deltaz * InvGradientVectorLength;
      {
      itk::Point<float, 3> BasePhysicalPoint;
      // std::cout<<__LINE__<<" :: "<<__FILE__<<std::endl;
      imapi->second->TransformIndexToPhysicalPoint(CurrentIndex, BasePhysicalPoint);
      // std::cout<<__LINE__<<" :: "<<__FILE__<<std::endl;
      for( float i = -ProfileExtent; i <= ProfileExtent; i += 1.0F )
        {
        // std::cout<<__LINE__<<" :: "<<__FILE__<< " :: i :: " <<i<<std::endl;

        itk::Point<float, 3> CurrProfilePoint = BasePhysicalPoint;

        CurrProfilePoint[0] = BasePhysicalPoint[0] + i * unitdeltax;
        CurrProfilePoint[1] = BasePhysicalPoint[1] + i * unitdeltay;
        CurrProfilePoint[2] = BasePhysicalPoint[2] + i * unitdeltaz;

        itk::ContinuousIndex<float, 3> ContinuousIndexProfilePoint;

        const bool inBuffer =
          imapi->second
          ->TransformPhysicalPointToContinuousIndex(CurrProfilePoint, ContinuousIndexProfilePoint);
        // From the imageinterpolators, Get the continuous index value.
        if( inBuffer )  //
          {
          // std::cout<<"MIN::"<<MapOfLocalMin[imapi->first] <<", "
          //         <<"MAX::"<<MapOfLocalMax[imapi->first] <<std::endl;
          *fi =  static_cast<neural_scalar_type>( Interpolator->
                                                  EvaluateAtContinuousIndex(ContinuousIndexProfilePoint) );

          /* linearly transform to zero one */
          // std::cout<<" (*fi - MapOfLocalMin[imapi->first] )::"
          //         <<(*fi - MapOfLocalMin[imapi->first] )<<std::endl;
          // std::cout<<"(MapOfLocalMax[imapi->first]-MapOfLocalMin[imapi->first] )::"
          //         <<(MapOfLocalMax[imapi->first]-MapOfLocalMin[imapi->first] )<<std::endl;
          // std::cout<<"*fi::"<<*fi;
          // std::cout<<" at::"<<ContinuousIndexProfilePoint<<std::endl;

          // normalization
          *fi = (*fi - MapOfLocalMin[imapi->first] ) / (MapOfLocalMax[imapi->first] - MapOfLocalMin[imapi->first] );

          // *fi = ( *fi -  MapOfLocalMax[imapi->first] )/ MapOfLocalMin[imapi->first] ;
          // std::cout<<__LINE__<<" :: "<<__FILE__<<":: "<< *fi<<std::endl;
          }
        else
          {
          std::stringstream ss;
          ss << "THIS IS HIGHLY UNLIKELY TO HAPPEN, AND NOT ALLOWED."
             << " The error is occuring at location " << CurrentIndex
             << " while accessing point " << ContinuousIndexProfilePoint << " with unit delta "
             << unitdeltax << " " << unitdeltay << " " << unitdeltaz << std::endl;
          ss << "Trying to access point :"
             << ContinuousIndexProfilePoint << " in image where physical space is from \n"
             << imapi->second->GetOrigin() << " to [";
          for( int q = 0; q < 3; q++ )
            {
            ss << imapi->second->GetOrigin()[q]
              + imapi->second->GetLargestPossibleRegion().GetSize()[q]
              * imapi->second->GetSpacing()[q] << ",";
            }
          ss << "]" << std::endl;
          std::cerr << ss.str();
          *fi = 0.0F;
          }

        // std::cout<<__LINE__<<" :: "<<__FILE__<<std::endl;
        fi++;
        // std::cout<<__LINE__<<" :: "<<__FILE__<<std::endl;
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

int AddSubjectInputVector(  DataSet * subjectSet,
                            NetConfiguration & ANNXMLObject,
                            const std::string registrationID,
                            const int inputVectorSize,
                            const int outputVectorSize,
                            const map<int, std::string>& MapOfROIOrder,
                            bool Apply) // last parameter: apply=false
{
  int VectorCreatedCounts = 0;

  std::cout << __LINE__ << " " << __FILE__ << std::endl;
  /** get necessary inputs
   * most of image identification will be done mapping from '-> Image' */

  /** image type to process */
  DataSet::StringVectorType ImageTypeList = ANNXMLObject.GetAtlasDataSet()->GetImageTypes();

  if( !ImageTypeList.size() )
    {
    itkGenericExceptionMacro(<< "No images types found. Cannot compute neural net output.");
    }

  /** read in images */
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  std::map<std::string, InternalImageType::Pointer> MapOfImages;
  for( DataSet::StringVectorType::const_iterator iIt = ImageTypeList.begin();
       iIt != ImageTypeList.end();
       ++iIt )
    {
    std::cout << *iIt << "::"
              << subjectSet->GetImageFilenameByType( *iIt) << std::endl;
    const std::string filename(  subjectSet->GetImageFilenameByType( *iIt) );
    MapOfImages[*iIt] = ReadImage<InternalImageType>( filename );
    }

  /** get registarion type */
  const RegistrationType *transform = subjectSet->GetRegistrationWithID(registrationID);
  std::string             DeformationFilename( transform->GetAttribute<StringValue>("AtlasToSubjRegistrationFilename") );

  /** read in Deformed probability map*/
  ProbabilityMapList * ROIList = ANNXMLObject.Get<ProbabilityMapList>("ProbabilityMapList");

  std::map<std::string, InternalImageType::Pointer> MapOfDeformedROIs;
  for( ProbabilityMapList::iterator ri = ROIList->begin();
       ri != ROIList->end();
       ++ri )
    {
    ProbabilityMapParser *ROI = dynamic_cast<ProbabilityMapParser *>(ri->second);
    MapOfDeformedROIs[ri->first] =
      ImageWarper<InternalImageType>( DeformationFilename,
                                      ROI->GetAttribute<StringValue>("Filename"),
                                      MapOfImages[*(ImageTypeList.begin() )] );
    }

  /** deform rho/phi/theta */
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  DataSet *                                         atlasDataSet = ANNXMLObject.GetAtlasDataSet();
  std::map<std::string, InternalImageType::Pointer> MapOfDeformedSpatialDescription;

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  MapOfDeformedSpatialDescription.insert(
    pair<std::string, InternalImageType::Pointer>(
      "rho", ImageWarper<InternalImageType>(
        DeformationFilename,
        atlasDataSet->GetSpatialLocationFilenameByType("rho"),
        MapOfImages[*(ImageTypeList.begin() )] ) ) );
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  MapOfDeformedSpatialDescription.insert(
    pair<std::string, InternalImageType::Pointer>(
      "phi", ImageWarper<InternalImageType>(
        DeformationFilename,
        atlasDataSet->GetSpatialLocationFilenameByType("phi"),
        MapOfImages[*(ImageTypeList.begin() )] ) ) );
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  MapOfDeformedSpatialDescription.insert(
    pair<std::string, InternalImageType::Pointer>(
      "theta", ImageWarper<InternalImageType>(
        DeformationFilename,
        atlasDataSet->GetSpatialLocationFilenameByType("theta"),
        MapOfImages[*(ImageTypeList.begin() )] ) ) );

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  NeuralParams *ANNmodel = ANNXMLObject.Get<NeuralParams>("NeuralNetParams");
  /** add ROIs vector */
  for( ProbabilityMapList::iterator pmi = ROIList->begin();
       pmi != ROIList->end();
       ++pmi )
    {
    std::cout << __LINE__ << "::" << __FILE__ << std::endl;
    if( !Apply )
      {
      std::cout << __LINE__ << "::" << __FILE__ << std::endl;
      /** ANNVectorStream */
      std::string   ANNVectorFilename = "InvalidANNVectors.txt";
      std::ofstream ANNVectorStream;
      // unsigned int  ANNTrainingVectorsCreated = 0;

      ANNVectorFilename = ANNmodel->GetAttribute<StringValue>("TrainingVectorFilename");
      ANNVectorFilename += "UnshuffledANN";

      std::string destination_dir =   itksys::SystemTools::GetFilenamePath(ANNVectorFilename);
      itksys::SystemTools::MakeDirectory( destination_dir.c_str() );

      ANNVectorStream.open(ANNVectorFilename.c_str(), std::ios::out | std::ios::binary | std::ios::app );

      std::cout << __LINE__ << "::" << __FILE__ << std::endl;
      if( !ANNVectorStream.good() )
        {
        itkGenericExceptionMacro(<< "Error: Could not open ANN vector file: "
                                 << ANNVectorFilename );
        }
      std::cout << __LINE__ << "::" << __FILE__ << std::endl;

      VectorCreatedCounts += AddROIVectorTrain( dynamic_cast<ProbabilityMapParser *>(pmi->second),
                                                (subjectSet),
                                                ANNXMLObject,
                                                MapOfROIOrder,
                                                MapOfImages,
                                                MapOfDeformedROIs,
                                                MapOfDeformedSpatialDescription,
                                                inputVectorSize,
                                                outputVectorSize,
                                                ANNVectorStream
                                                );
      ANNVectorStream.close();
      }
    else
      {
      std::cout << __LINE__ << "::" << __FILE__ << std::endl;
      /** get the trained model file name */
      std::string ANNModelFilename = ANNmodel->GetAttribute<StringValue>("TrainingModelFilename");
      ANNModelFilename += "ANN";
      ANNParams *annParams = ANNXMLObject.Get<ANNParams>("ANNParams");
      int        TimePoint      = annParams->GetAttribute<IntValue>("Iterations");

      /** check if optimally trained point is given */
      std::string ANNVectorHeaderFilename = ANNmodel->GetAttribute<StringValue>("TrainingVectorFilename");
      ANNVectorHeaderFilename += "ANN.hdr";

      std::ifstream headerStream;
      headerStream.open( ANNVectorHeaderFilename.c_str(), std::ios::in | std::ios::binary);

      std::cout << __LINE__ << "::" << __FILE__ << std::endl;
      char currentline[MAX_LINE_SIZE];

      std::cout << __LINE__ << "::" << __FILE__ << std::endl;
      if( !headerStream.is_open() )
        {
        std::cout << "Error: Fail to open file of "
                  << ANNVectorHeaderFilename
                  << std::endl;
        std::cout << "Use given iteration number" << std::endl;
        }
      else
        {
        while( headerStream.getline(currentline, MAX_LINE_SIZE - 1) )
          {
          std::cout << __LINE__ << "::" << __FILE__ << std::endl;
          std::string        temp;
          std::istringstream iss(currentline, std::istringstream::in);
          iss >> temp;
          if( temp == ANNModelFilename )
            {
            iss >> TimePoint;
            std::cout << "================================================="
                      << std::endl
                      << "** Using optimally trained point "
                      << TimePoint << std::endl;
            std::cout << "================================================="
                      << std::endl;
            }
          }
        }

      char temp[10];
      sprintf( temp, "%09d", TimePoint );
      ANNModelFilename += temp;

      std::cout << __LINE__ << "::" << __FILE__ << std::endl;
      AddROIVectorApply( dynamic_cast<ProbabilityMapParser *>(pmi->second),
                         (subjectSet),
                         ANNXMLObject,
                         MapOfROIOrder,
                         MapOfImages,
                         MapOfDeformedROIs,
                         MapOfDeformedSpatialDescription,
                         inputVectorSize,
                         outputVectorSize,
                         ANNModelFilename
                         );
      std::cout << __LINE__ << "::" << __FILE__ << std::endl;
      }
    }

  return VectorCreatedCounts;
}

int AddROIVectorTrain( ProbabilityMapParser * currentROI,
                       DataSet *subjectSet,
                       NetConfiguration & ANNXMLObject,
                       std::map<int, std::string>  MapOfROIOrder,
                       std::map<std::string, InternalImageType::Pointer> & MapOfImages,
                       std::map<std::string, InternalImageType::Pointer> & MapOfDeformedROIs,
                       std::map<std::string, InternalImageType::Pointer> & MapOfDeformedSpatialDescription,
                       const int inputVectorSize,
                       const int outputVectorSize,
                       std::ofstream & ANNVectorStream
                       )
{
  /** create vectors for current ROI */

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  std::string currentROIName( currentROI->GetAttribute<StringValue>("StructureID") );

  InternalImageType::Pointer knownOutputVolume;

  std::string knownBinaryVolumeName( subjectSet->GetMaskFilenameByType( currentROIName) );

  std::cout << __LINE__ << "::" << __FILE__ << "::knownBinaryVolumeName::" << knownBinaryVolumeName << std::endl;

  knownOutputVolume
    = itkUtil::ReadImage<InternalImageType>( knownBinaryVolumeName);
  std::cout << __LINE__ << "::" << __FILE__ << "::knownOutputVolume::" << knownOutputVolume << std::endl;

  /** get the local mean and std */
  std::map<std::string, neural_scalar_type> MapOfLocalMax;
  std::map<std::string, neural_scalar_type> MapOfLocalMin;

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  DataSet::StringVectorType ImageTypeList = ANNXMLObject.GetAtlasDataSet()->GetImageTypes();
  if( !ImageTypeList.size() )
    {
    itkGenericExceptionMacro(<< "No images types found. Cannot compute neural net output.");
    }
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  // create labelmap from probability

  typedef itk::Image<unsigned char, 3> BinaryImageType;

  typedef itk::BinaryThresholdImageFilter<InternalImageType,
                                          BinaryImageType> BinaryThresholdType;

  BinaryThresholdType::Pointer thresholder = BinaryThresholdType::New();

  thresholder->SetInput( MapOfDeformedROIs[currentROIName] );
  thresholder->SetLowerThreshold( 0.001F );
  thresholder->SetInsideValue(1);
  thresholder->Update();

  itkUtil::WriteImage<BinaryImageType>( thresholder->GetOutput(),
                                        "./DEBUG_BinaryImage" + currentROIName + ".nii.gz");

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  for( DataSet::StringVectorType::const_iterator ImageList = ImageTypeList.begin();
       ImageList != ImageTypeList.end();
       ++ImageList )
    {
    /**[ local normalization ]
    * The local normalization process is done with mean and standard deviation
    * of intensity for the region. The region is identified with deformed
    * probability map where it has a non-zero probability */

    // calculate mean and std in the ROI defined by binary mask
    itkUtil::WriteImage<InternalImageType>( MapOfImages[*ImageList],
                                            "./DEBUG_MapOfImage.nii.gz");
    typedef itk::LabelStatisticsImageFilter<InternalImageType,
                                            BinaryImageType> StatCalculatorType;

    StatCalculatorType::Pointer statCalculator = StatCalculatorType::New();
    statCalculator->SetInput( MapOfImages[*ImageList] );
    statCalculator->SetLabelInput( thresholder->GetOutput() );
    statCalculator->Update();

    // get the mean and std
    MapOfLocalMax[*ImageList] = statCalculator->GetMaximum( thresholder->GetInsideValue() );
    MapOfLocalMin[*ImageList] = statCalculator->GetMinimum( thresholder->GetInsideValue() );
    std::cout << __LINE__ << "::" << __FILE__ << "::Max::" << MapOfLocalMax[*ImageList] << std::endl;
    std::cout << __LINE__ << "::" << __FILE__ << "::Min::" << MapOfLocalMin[*ImageList] << std::endl;

    /*
    itk::ImageRegionIterator<InternalImageType> ti( MapOfImages[*ImageList] ,
                                                    MapOfImages[*ImageList]->GetLargestPossibleRegion());
    ti.GoToBegin();
    float max,min;
    max= MapOfImages[*ImageList]->GetPixel( ti.GetIndex() );
    min= MapOfImages[*ImageList]->GetPixel( ti.GetIndex() );
    while(!ti.IsAtEnd() )
    {
      InternalImageType::IndexType temp = ti.GetIndex();
      if( MapOfDeformedROIs[ currentROIName ]->GetPixel( temp ) > 0.001F)
      {
          if( MapOfImages[*ImageList]->GetPixel( temp ) > max)
          {
            max = MapOfImages[*ImageList]->GetPixel( temp ) ;
          }
          else
          {
            min=MapOfImages[*ImageList]->GetPixel( temp ) ;
          }
      }
      ++ti;
    }
    std::cout<<__LINE__<<"::"<<__FILE__<<"::Max1::"<<max<<std::endl;
    std::cout<<__LINE__<<"::"<<__FILE__<<"::Min2::"<<min<<std::endl;
*/
    }

  /** get gradient profile size */
  NeuralParams *ANNmodel = ANNXMLObject.Get<NeuralParams>("NeuralNetParams");
  const int     GradientProfileSize = ANNmodel->GetAttribute<IntValue>("GradientProfileSize");

  /** count region */
  itk::Index<3> min, max;
  DefineBoundingBox(MapOfDeformedROIs[currentROIName], min, max);

  InternalImageType::RegionType BBRegion;
  itk::Size<3>                  BBRange;
  itk::Index<3>                 BBIndex;

  unsigned int long range[3] = { max[0] - min[0], max[1] - min[1], max[2] - min[2]};
  long int          index[3] = { min[0], min[1], min[2]};

  BBRange.SetSize(range);
  BBRegion.SetSize(BBRange);
  BBIndex.SetIndex(index);
  BBRegion.SetIndex(BBIndex);
  /** count non zero probablility region */
  unsigned int NonZeroProbMapPoints = 0;

  itk::ImageRegionIterator<InternalImageType> bbri(MapOfDeformedROIs[currentROIName], BBRegion);
  while( !bbri.IsAtEnd() )
    {
    const InternalImageType::IndexType CurrentIndex = bbri.GetIndex();
    InternalImageType::PixelType       current_value =
      MapOfDeformedROIs[currentROIName]->GetPixel(CurrentIndex);
    if( ( current_value < ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) ) &&
        ( current_value > ( 0.0F + PERCENT_MIN_MAX_TOLERANCE ) ) )
      {
      NonZeroProbMapPoints++;
      }
    ++bbri;
    }

  if( NonZeroProbMapPoints == 0 )
    {
    itkGenericExceptionMacro(<< "ERROR:  Could not find any valid probable points for deformed ROI");
    }

  /** gradient direction image */
  itk::GradientImageFilter<InternalImageType,
                           float,
                           float>::Pointer GradientFilter =
    itk::GradientImageFilter<InternalImageType, float, float>::New();

  GradientFilter->SetInput(MapOfDeformedROIs[currentROIName]);
  GradientFilter->Update();
  itk::Image<itk::CovariantVector<float,
                                  3>, 3>::Pointer gradient = GradientFilter->GetOutput();

  /** input/output vector allocator */
  std::vector<neural_scalar_type> inputvector(inputVectorSize);
  std::vector<neural_scalar_type> outputvector(outputVectorSize);

  const unsigned int HardMaxNumberSamplingPoints = 500000;

  int SamplesLeftToFind =
    ( HardMaxNumberSamplingPoints < NonZeroProbMapPoints ) ?
    HardMaxNumberSamplingPoints : NonZeroProbMapPoints;
  float *write_buffer =
    new float[SamplesLeftToFind * (inputVectorSize + outputVectorSize + 1)];

  int value_index = 0;
  int currentBufferCounter = 0;

  /** random walk iterator
   *  set the number of sample as a whole region and then
   *  count how many voxels have been inserted as a vector. */
  itk::ImageRandomConstIteratorWithIndex<InternalImageType>
  randomWalkRegionIterator( MapOfDeformedROIs[currentROIName], BBRegion );

  randomWalkRegionIterator.GoToBegin();
  randomWalkRegionIterator.SetNumberOfSamples( BBRegion.GetNumberOfPixels() );

  // std::cout<<__LINE__<<__FILE__<<std::endl;
  while(  (!randomWalkRegionIterator.IsAtEnd() ) &&
          (SamplesLeftToFind) )
    {
    const InternalImageType::IndexType currentIndex = randomWalkRegionIterator.GetIndex();

    InternalImageType::PixelType current_value =
      static_cast<InternalImageType::PixelType>(
        MapOfDeformedROIs[currentROIName]->GetPixel( currentIndex ) );

    neural_scalar_type guard(AUTOSEG_VEC_SENTINEL);

    if( ( current_value < ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) ) &&
        ( current_value > ( 0.0F + PERCENT_MIN_MAX_TOLERANCE ) ) )
      {
      AddInputVector( inputvector,
                      gradient,
                      MapOfDeformedSpatialDescription["rho"],
                      MapOfDeformedSpatialDescription["phi"],
                      MapOfDeformedSpatialDescription["theta"],
                      MapOfImages,
                      MapOfDeformedROIs,
                      MapOfLocalMax,
                      MapOfLocalMin,
                      MapOfROIOrder,
                      currentIndex,
                      GradientProfileSize);

      /** add output vector */
      std::cout << "Index" << currentIndex << ", ";
      std::vector<neural_scalar_type>::iterator oi = outputvector.begin();
      for( unsigned int i = 0; i < outputvector.size(); i++ )
        {
        if( MapOfROIOrder.find(i)->second == currentROIName )
          {
          if( knownOutputVolume->GetPixel(currentIndex) > 0.5F )
            {
            *oi = 1.0F;
            }
          else
            {
            *oi = 0.0F;
            }
          }
        else
          {
          *oi = 0.0F;
          }
        ++oi;
        }
      // writing out outputVectorFirst
      for( std::vector<neural_scalar_type>::const_iterator ofi = outputvector.begin();
           ofi != outputvector.end();
           ++ofi )
        {
        std::cout << *ofi << ", ";
        write_buffer[value_index] = *ofi;
        value_index++;
        }
      for( std::vector<neural_scalar_type>::const_iterator ifi = inputvector.begin();
           ifi != inputvector.end();
           ++ifi )
        {
        std::cout << *ifi << ", ";
        write_buffer[value_index] = *ifi;
        value_index++;
        }
      // std::cout<<endl;
      /** add guard */
      write_buffer[value_index] = guard;
      value_index++;

      currentBufferCounter++;
      SamplesLeftToFind--;
      }
    ++randomWalkRegionIterator;
    }

  ANNVectorStream.write( (const char *)( write_buffer ),
                         sizeof( float ) * (inputVectorSize + outputVectorSize + 1) * currentBufferCounter );
  delete[] write_buffer;
  return currentBufferCounter;
}

void AddROIVectorApply( ProbabilityMapParser * currentROI,
                        DataSet *subjectSet,
                        NetConfiguration & ANNXMLObject,
                        std::map<int, std::string> MapOfROIOrder,
                        std::map<std::string, InternalImageType::Pointer>& MapOfImages,
                        std::map<std::string, InternalImageType::Pointer>& MapOfDeformedROIs,
                        std::map<std::string, InternalImageType::Pointer>& MapOfDeformedSpatialDescription,
                        const int inputVectorSize,
                        const int outputVectorSize,
                        const std::string  ANNModelFilename
                        )
{
  /** create vectors for current ROI */

  std::string currentROIName( currentROI->GetAttribute<StringValue>("StructureID") );

  int currentROINumber = -1;

  for( int i = 0; i < outputVectorSize; i++ )
    {
    if( MapOfROIOrder.find(i)->second == currentROIName )
      {
      currentROINumber = i;
      }
    }
  if( currentROINumber == -1 )
    {
    itkGenericExceptionMacro(<< "No valid ROINumber found" );
    }

  /** get the local mean and std */
  std::map<std::string, neural_scalar_type> MapOfLocalMax;
  std::map<std::string, neural_scalar_type> MapOfLocalMin;

  DataSet::StringVectorType ImageTypeList = ANNXMLObject.GetAtlasDataSet()->GetImageTypes();
  if( !ImageTypeList.size() )
    {
    itkGenericExceptionMacro(<< "No images types found. Cannot compute neural net output.");
    }
  // create labelmap from probability

  typedef itk::Image<unsigned char, 3> BinaryImageType;

  typedef itk::BinaryThresholdImageFilter<InternalImageType,
                                          BinaryImageType> BinaryThresholdType;

  BinaryThresholdType::Pointer thresholder = BinaryThresholdType::New();

  thresholder->SetInput( MapOfDeformedROIs[currentROIName] );
  thresholder->SetLowerThreshold( 0.0001F );
  thresholder->SetInsideValue(1);
  thresholder->Update();

  itkUtil::WriteImage<BinaryImageType>( thresholder->GetOutput(),
                                        "./DEBUG_BinaryImage" + currentROIName + ".nii.gz");
  for( DataSet::StringVectorType::const_iterator ImageList = ImageTypeList.begin();
       ImageList != ImageTypeList.end();
       ++ImageList )
    {
    /**[ local normalization ]
    * The local normalization process is done with mean and standard deviation
    * of intensity for the region. The region is identified with deformed
    * probability map where it has a non-zero probability */

    // calculate mean and std in the ROI defined by binary mask
    itkUtil::WriteImage<InternalImageType>( MapOfImages[*ImageList],
                                            "./DEBUG_MapOfImage.nii.gz");
    typedef itk::LabelStatisticsImageFilter<InternalImageType,
                                            BinaryImageType> StatCalculatorType;

    StatCalculatorType::Pointer statCalculator = StatCalculatorType::New();

    statCalculator->SetInput( MapOfImages[*ImageList] );
    statCalculator->SetLabelInput( thresholder->GetOutput() );
    statCalculator->Update();

    // get the mean and std
    MapOfLocalMax[*ImageList] = statCalculator->GetMaximum( thresholder->GetInsideValue() );
    MapOfLocalMin[*ImageList] = statCalculator->GetMinimum( thresholder->GetInsideValue() );
    std::cout << __LINE__ << "::" << __FILE__ << "::Max::" << MapOfLocalMax[*ImageList] << std::endl;
    std::cout << __LINE__ << "::" << __FILE__ << "::Min::" << MapOfLocalMin[*ImageList] << std::endl;
    /*
     itk::ImageRegionIterator<InternalImageType> ti( MapOfImages[*ImageList] ,
                                                     MapOfImages[*ImageList]->GetLargestPossibleRegion());
     ti.GoToBegin();
     float max,min;
     max= MapOfImages[*ImageList]->GetPixel( ti.GetIndex() );
     min= MapOfImages[*ImageList]->GetPixel( ti.GetIndex() );
     while(!ti.IsAtEnd() )
     {
       InternalImageType::IndexType temp = ti.GetIndex();
       if( MapOfDeformedROIs[ currentROIName ]->GetPixel( temp ) > 0.001F)
       {
           if( MapOfImages[*ImageList]->GetPixel( temp ) > max)
           {
             max = MapOfImages[*ImageList]->GetPixel( temp ) ;
           }
           else
           {
             min=MapOfImages[*ImageList]->GetPixel( temp ) ;
           }
       }
       ++ti;
     }
     std::cout<<__LINE__<<"::"<<__FILE__<<"::Max1::"<<max<<std::endl;
     std::cout<<__LINE__<<"::"<<__FILE__<<"::Min2::"<<min<<std::endl;
   */
    }

  /** get gradient profile size */
  NeuralParams *ANNmodel = ANNXMLObject.Get<NeuralParams>("NeuralNetParams");
  const int     GradientProfileSize = ANNmodel->GetAttribute<IntValue>("GradientProfileSize");

  /** count region */
  itk::Index<3> min, max;
  DefineBoundingBox(MapOfDeformedROIs[currentROIName], min, max);

  InternalImageType::RegionType BBRegion;
  itk::Size<3>                  BBRange;
  itk::Index<3>                 BBIndex;

  unsigned int long range[3] = { max[0] - min[0], max[1] - min[1], max[2] - min[2]};
  long int          index[3] = { min[0], min[1], min[2]};

  BBRange.SetSize(range);
  BBRegion.SetSize(BBRange);
  BBIndex.SetIndex(index);
  BBRegion.SetIndex(BBIndex);
  /** count non zero probablility region */
  unsigned int NonZeroProbMapPoints = 0;

  itk::ImageRegionIterator<InternalImageType> bbri(MapOfDeformedROIs[currentROIName], BBRegion);
  bbri.GoToBegin();
  while( !bbri.IsAtEnd() )
    {
    const InternalImageType::IndexType CurrentIndex = bbri.GetIndex();
    InternalImageType::PixelType       current_value =
      MapOfDeformedROIs[currentROIName]->GetPixel(CurrentIndex);
    if( ( current_value < ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) ) &&
        ( current_value > ( 0.0F + PERCENT_MIN_MAX_TOLERANCE ) ) )
      {
      NonZeroProbMapPoints++;
      }
    ++bbri;
    }

  if( NonZeroProbMapPoints == 0 )
    {
    itkGenericExceptionMacro(<< "ERROR:  Could not find any valid probable points for deformed ROI");
    }

  /** gradient direction image */
  itk::GradientImageFilter<InternalImageType,
                           float,
                           float>::Pointer GradientFilter =
    itk::GradientImageFilter<InternalImageType, float, float>::New();

  GradientFilter->SetInput(MapOfDeformedROIs[currentROIName]);
  GradientFilter->Update();
  itk::Image<itk::CovariantVector<float,
                                  3>, 3>::Pointer gradient = GradientFilter->GetOutput();

  /** predicted value storage */
  typedef std::list<neural_scalar_type>           predictedValueListType;
  typedef std::list<InternalImageType::IndexType> predictedIndexListType;

  predictedValueListType predictedValueList;
  predictedIndexListType predictedIndexList;

  /** read in trained model */
  neural_net_type * ANNModel = new neural_net_type();

  std::cout << __LINE__ << "::" << __FILE__ << "::" << ANNModelFilename << std::endl;
  if( !itksys::SystemTools::FileExists( ANNModelFilename.c_str() ) )
    {
    std::cout << __LINE__ << "::" << __FILE__ << std::endl
              << " Model file dose not exist :: "
              << ANNModelFilename << std::endl;
    }
  ANNModel->load( ANNModelFilename.c_str() );

  /** input/output vector allocator */
  std::vector<neural_scalar_type> inputvector(inputVectorSize);
  std::vector<neural_scalar_type> outputvector(outputVectorSize);
  bbri.GoToBegin();

  /** get the limitaion */
  while(  !bbri.IsAtEnd() )
    {
    neural_scalar_type                 tempPrediction;
    const InternalImageType::IndexType currentIndex = bbri.GetIndex();

    InternalImageType::PixelType current_value = static_cast<InternalImageType::PixelType>(
        MapOfDeformedROIs[currentROIName]->GetPixel( currentIndex ) );
    // std::cout<<__LINE__<<" :: "<<__FILE__<<" :: current_value :: " << current_value<<std::endl;

    if( current_value <  (0.0F + PERCENT_MIN_MAX_TOLERANCE) )
      {
      tempPrediction = 0.0F;
      // std::cout<<"::tempPrediction::"<<tempPrediction<<std::endl;
      }
    else if( current_value < ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) )
      {
      AddInputVector( inputvector,
                      gradient,
                      MapOfDeformedSpatialDescription["rho"],
                      MapOfDeformedSpatialDescription["phi"],
                      MapOfDeformedSpatialDescription["theta"],
                      MapOfImages,
                      MapOfDeformedROIs,
                      MapOfLocalMax,
                      MapOfLocalMin,
                      MapOfROIOrder,
                      currentIndex,
                      GradientProfileSize);
      /* predict based on the input vector */

      neural_vector_type inputInOpenCV = cvCreateMat( 1, inputVectorSize, CV_32FC1);
      neural_vector_type outputInOpenCV = cvCreateMat( 1, outputVectorSize, CV_32FC1);

      neural_scalar_type *tempInput = new neural_scalar_type[inputVectorSize];

      std::vector<float>::const_iterator fi = inputvector.begin();
      for( int j = 0; j < inputVectorSize; j++, fi++ )
        {
        // std::cout<<*fi<<", ";
        tempInput[j] = *fi;
        }
      cvInitMatHeader( inputInOpenCV, 1, inputVectorSize, CV_32FC1, tempInput);

      ANNModel->predict( inputInOpenCV, outputInOpenCV );

      tempPrediction = CV_MAT_ELEM( *outputInOpenCV, neural_scalar_type, 0, currentROINumber);

      tempPrediction *= HUNDRED_PERCENT_VALUE;
      if( tempPrediction > HUNDRED_PERCENT_VALUE )
        {
        tempPrediction = HUNDRED_PERCENT_VALUE;
        }
      else if( tempPrediction < 0.0F )
        {
        tempPrediction = 0.0F;
        }

      // std::cout<<"::tempPrediction::"<<tempPrediction<<std::endl;
      }
    else
      {
      tempPrediction = HUNDRED_PERCENT_VALUE;
      // std::cout<<"::tempPrediction::"<<tempPrediction<<std::endl;
      }

    predictedValueList.push_back( tempPrediction);
    predictedIndexList.push_back( currentIndex);

    ++bbri;
    }

  /** predicted output image */
  InternalImageType::Pointer ANNOutput = InternalImageType::New();

  ANNOutput->CopyInformation( MapOfImages.begin()->second );
  ANNOutput->SetRegions( MapOfImages.begin()->second->GetLargestPossibleRegion() );
  ANNOutput->Allocate();
  ANNOutput->FillBuffer( 0.0F );

  predictedIndexListType::const_iterator pIndexi = predictedIndexList.begin();
  for( predictedValueListType::const_iterator pi = predictedValueList.begin();
       pi != predictedValueList.end();
       pi++, pIndexi++ )
    {
    ANNOutput->SetPixel( *pIndexi, *pi );
    }

  /** write out Predicted Output */
  std::string outputDir = subjectSet->GetAttribute<StringValue>("OutputDir");
  itksys::SystemTools::MakeDirectory( outputDir.c_str() );

  const std::string subjectID(subjectSet->GetAttribute<StringValue>("Name") );

  std::string outputVolumeFilename = outputDir + "/ANNContinuousPrediction"
    + currentROIName
    + subjectID
    + ".nii.gz";

  itkUtil::WriteImage<InternalImageType>( ANNOutput, outputVolumeFilename );

  /** write out Median Filtered Predicted Output */
  typedef itk::MedianImageFilter<InternalImageType,
                                 InternalImageType> MedianFilterType;

  MedianFilterType::Pointer medianFilter = MedianFilterType::New();

  InternalImageType::SizeType RadiusType;

  RadiusType.Fill(1);

  medianFilter->SetRadius(RadiusType);
  medianFilter->SetInput( ANNOutput );
  medianFilter->Update();

  std::cout << __LINE__ << "::" << __FILE__
            << "::CREATE OUTPUT :: "
            <<
    outputDir
    + "/ANNContinuousPredictionMedianFiltered_"
    + currentROIName
    + subjectID
    + ".nii.gz" << std::endl;
  itkUtil::WriteImage<InternalImageType>( medianFilter->GetOutput(),
                                          outputDir
                                          + "/ANNContinuousPredictionMedianFiltered_"
                                          + currentROIName
                                          + subjectID
                                          + ".nii.gz"
                                          );
}

void AddInputVector
  (std::vector<neural_scalar_type> & inputvector,
  itk::Image<itk::CovariantVector<float, 3>, 3>::Pointer ProbMapGradient,
  InternalImageType::Pointer RhoMapImage,
  InternalImageType::Pointer PhiMapImage,
  InternalImageType::Pointer ThetaMapImage,
  const std::map<std::string, InternalImageType::Pointer> & MapOfImages,
  const std::map<std::string, InternalImageType::Pointer> & MapOfDeformedROIs,
  std::map<std::string, neural_scalar_type>& MapOfLocalMax,
  std::map<std::string, neural_scalar_type>& MapOfLocalMin,
  std::map<int, std::string>& MapOfROIOrder,
  const InternalImageType::IndexType & CurrentIndex,
  const int GradientProfileSize)
{
  // std::cout<<"Add input Vector...."<<std::endl;
  std::vector<neural_scalar_type>::iterator fi = inputvector.begin();
  // std::cout<<"EMbrased part start here.....:"<<std::endl;
    {
    // std::cout<< __FILE__ << " " << __LINE__ << std::endl;
    *fi = 1.0F;   // Hardcode to 1.0 just to eliminate it's effect
    // static_cast<float>(ProbMapImage->GetPixel(CurrentIndex));
    // *fi=static_cast<float>(ProbMapImage->GetPixel(CurrentIndex));
    // std::cout<<"\n"<<*fi<<",";
    // fi++;
    /*========================================================================//
      * Add Sphericcal Coord System
      * =========================================================================*/
    *fi = static_cast<neural_scalar_type>( RhoMapImage->GetPixel(CurrentIndex) );
    ++fi;
    *fi = static_cast<neural_scalar_type>( PhiMapImage->GetPixel(CurrentIndex) );
    ++fi;
    *fi = static_cast<neural_scalar_type>( ThetaMapImage->GetPixel(CurrentIndex) );
    ++fi;
    /*========================================================================//
      * Add Probability Maps....
      * =========================================================================*/
    // std::cout<< __FILE__ << " " << __LINE__ << std::endl;
    for( unsigned int i = 0; i < MapOfROIOrder.size(); i++ )
      {
      std::string                roiName( MapOfROIOrder.find(i)->second);
      InternalImageType::Pointer temp = MapOfDeformedROIs.find( roiName )->second;
      float                      probabilityValue =
        ( MapOfDeformedROIs.find( roiName )->second)->GetPixel(CurrentIndex);

      if( probabilityValue > 0.0 )
        {
        *fi = 1.0F;
        }
      else
        {
        *fi = 0.0F;
        }
      ++fi;
      }
    // std::cout<< __FILE__ << " " << __LINE__ << std::endl;
    }

  /*========================================================================//
    * Find Gradient Value....
    * =========================================================================*/
  // std::cout<< __FILE__ << " " << __LINE__ << std::endl;

  FillGradProfile(fi,                           // Index + rho phi theta
                  MapOfImages,
                  MapOfLocalMax,
                  MapOfLocalMin,
                  ProbMapGradient,
                  CurrentIndex,
                  GradientProfileSize/*,
                                         * IrisSize--> We are not using iris size
                                         * any more.  */);
  // std::cout<< __FILE__ << " " << __LINE__ << std::endl;
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

TDeformationField::Pointer
RegisterLandmarksToDeformationField(const std::string & InputImageFilename,
                                    const std::string & InputAtlasFilename,
                                    const std::string & TemplateAtlasFilename,
                                    const std::string & TemplateImageFilename)
{
  return itkUtil::RegisterLandmarksToDeformationField
         <TDeformationField, InternalImageType>(InputImageFilename,
                                                InputAtlasFilename,
                                                TemplateAtlasFilename,
                                                TemplateImageFilename);
}

#endif
