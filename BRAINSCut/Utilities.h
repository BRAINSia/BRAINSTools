/**
  * \mainpage Artificial Neural Network Image Segmentation Documentation
  * \section intro_sec Introduction
  * Why did we do this?
  */
#ifndef __Utilities_h__
#define __Utilities_h__

/** The following definitions are temporaries used to help identify code that
  * must be refactored in order to remove all
  * dependancies of probability images being written as un-signed characters's
  * */

static const double HUNDRED_PERCENT_VALUE = 1.0F;
static const double MINIMUM_VALUE = 0.0F;
static const double PERCENT_MIN_MAX_TOLERANCE = 0.01F;
// static const float
// THRESHOLD_CUTOFF=HUNDRED_PERCENT_VALUE-PERCENT_MIN_MAX_TOLERANCE;
// static const float
// INV_HUNDRED_PERCENT_VALUE=1.0F/static_cast<float>(HUNDRED_PERCENT_VALUE);

#include <stdint.h>
// code consolidation for reading and orienting and writing and orienting the
// image.
#include "itkIO.h"
#include "itksys/SystemTools.hxx"
#include "itksys/Process.h"

#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

// ITK
#include "itkAddImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkBatchSupervisedTrainingFunction.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkDifferenceOfGaussiansGradientImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkElasticBodySplineKernelTransform.h"
#include "itkGradientImageFilter.h"
#include "itkHausdorffDistanceImageFilter.h"
#include "itkImage.h"
#include "itkImageIOFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegistrationMethod.h"
#include "itkImportImageFilter.h"
#include "itkIOCommon.h"
#include "itkIterativeSupervisedTrainingFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkListSample.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkMetaImageIO.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkNumericTraits.h"
#include "itkOneHiddenLayerBackPropagationNeuralNetwork.h"
#include "itkOrientImageFilter.h"
#include "itkPoint.h"
#include "itkPointSet.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkReinitializeLevelSetImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSimilarityIndexImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "itkThresholdImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkVector.h"
#include "itkWarpImageFilter.h"
#include "itkWeightSetBase.h"

#include "GenericTransformImage.h"
#include "itkLabelStatisticsImageFilter.h"

// ****************OpenCV****************//
// ****************OpenCV****************//
#define USE_OPENCV 1##NEED FOR SVM
#include "ml.h"
#include "cxcore.h"
#include "cv.h"
typedef float   neural_scalar_type;
typedef CvMat * neural_vector_type;

typedef struct
  {
  neural_vector_type inputVector;
  neural_vector_type outputVector;
  int numberOfVector;
  } neural_data_set_type;

typedef CvANN_MLP_Revision neural_net_type;
// typedef CvANN_MLP neural_net_type;
// **************************************//

// ******************SVNM****************//
// ******************SVNM****************//
#include "Svm.h"
#include "NetConfiguration.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_c_vector.h"
#include <map>
#include <list>
#include <fstream>
#include <float.h>
#include <iostream>
#include <ctime>
#include <exception>
#include <sstream>
#include <itksys/SystemTools.hxx>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <stdarg.h>
#include <signal.h>
#include <stddef.h>
#include <ctype.h>
#include <DataSet.h>

#if defined( _WIN32 ) && ( defined( _MSC_VER ) || defined( __BORLANDC__ ) )
#include <sys\timeb.h>
#include <dos.h>
#define _unlink unlink
#else
#include <unistd.h>
#include <sys/timeb.h>
#endif

#if !defined( MAX_LINE_SIZE )
#define MAX_LINE_SIZE 4096
#endif

#define AUTOSEG_VEC_SENTINEL 1234567.0
#if defined( USE_BRAINS2_TRANSFORMS )
#include "b2Affine_rw.h"
#endif // USE_BRAINS2_TRANSFORMS

//
// BRAINS HELPER
//
#include "BRAINSFitHelper.h"

typedef itk::Image<float, 3> InputImageType;
typedef itk::Image<float, 3> InternalImageType;
typedef itk::Image<itk::Vector<float,
                               3>, 3>                     TDeformationField;
typedef itk::Image<itk::CovariantVector<float,
                                        3>, 3>                     GradImageType;
typedef itk::RecursiveGaussianImageFilter<itk::Image<float,
                                                     3>, itk::Image<float, 3> > g_GaussianFilterType;

#if 0
const unsigned int WindowRadius = 5;
typedef itk::Function::HammingWindowFunction<WindowRadius, float,
                                             float> WindowFunctionType;

typedef itk::ConstantBoundaryCondition<InternalImageType>
  BoundaryConditionType;
typedef itk::WindowedSincInterpolateImageFunction<
    InternalImageType,
    WindowRadius,
    WindowFunctionType,
    BoundaryConditionType,
    float>   ImageLinearInterpolatorType;

typedef itk::ConstantBoundaryCondition<InputImageType> RealBoundaryConditionType;
typedef itk::WindowedSincInterpolateImageFunction<
    InputImageType,
    WindowRadius,
    WindowFunctionType,
    RealBoundaryConditionType,
    float>  ImageLinearInterpolatorRealType;
#else
typedef itk::LinearInterpolateImageFunction<InternalImageType,
                                            float> ImageLinearInterpolatorType;
typedef itk::LinearInterpolateImageFunction<InputImageType,
                                            float> ImageLinearInterpolatorRealType;
#endif

extern bool CHECK_CORONAL(InternalImageType::DirectionType Dir);

/**
  * \ingroup Util
  * \short
  *   This function smooths and image and writes the image to disk.
  * \param InputImageFilename Input Image Filename
  * \param currentGaussianSigma Gauassian Smoothing Value
  * \param OutputImageFilename Output Image Filename
  * \return Void
  */
extern void SmoothImageWrite(const std::string & InputImageFilename, const float currentGaussianSigma,
                             const std::string & OutputImageFilename);

/**
  * \ingroup Util
  * \short
  *   This function reads a deformation field.
  * \param DeformationField Deformation Field Poitner
  * \param DeformationFilename Input Deformation Filename
  * \return Void
  */
extern void ReadDeformationField(TDeformationField::Pointer & DeformationField,
                                 const std::string & DeformationFilename);

/**
  * \ingroup Util
  * \short
  *  This function smooths a single image.
  * \param ImageToSmooth Input Image Pointer
  * \param currentGaussianSigma Gaussian Smoothing Value
  * \param SmoothedImage Output Image Pointer
  * \return Void
  */
// extern void SmoothSingleImage(InternalImageType::Pointer
// &ImageToSmooth,
//     const float currentGaussianSigma,
//     InternalImageType::Pointer &SmoothedImage);
template <class SmootherImageType>
typename SmootherImageType::Pointer
SmoothSingleImage(typename SmootherImageType::Pointer ImageToSmooth,
                  const float currentGaussianSigma)
{
  if( currentGaussianSigma <= 0.0F )
    {
    return ImageToSmooth;
    }
  try
    {
    std::cout << " Smoothing Guassian with of " << currentGaussianSigma
              << std::endl;
    typename itk::RescaleIntensityImageFilter<SmootherImageType,
                                              InputImageType>::Pointer
    ProbMapToRealCast =
      itk::RescaleIntensityImageFilter<SmootherImageType, InputImageType>::New();
    ProbMapToRealCast->SetOutputMinimum(0.0F);
    ProbMapToRealCast->SetOutputMaximum(1.0F);
    ProbMapToRealCast->SetInput(ImageToSmooth);
    ProbMapToRealCast->ReleaseDataFlagOn();
    typename g_GaussianFilterType::Pointer g_gaussianFilterX =
      g_GaussianFilterType::New();
    g_gaussianFilterX->SetNormalizeAcrossScale(true);
    g_gaussianFilterX->SetSigma(currentGaussianSigma);
    g_gaussianFilterX->SetDirection(0);
    g_gaussianFilterX->SetInput( ProbMapToRealCast->GetOutput() );
    g_gaussianFilterX->ReleaseDataFlagOn();
    typename g_GaussianFilterType::Pointer g_gaussianFilterY =
      g_GaussianFilterType::New();
    g_gaussianFilterY->SetNormalizeAcrossScale(true);
    g_gaussianFilterY->SetSigma(currentGaussianSigma);
    g_gaussianFilterY->SetDirection(1);
    g_gaussianFilterY->SetInput( g_gaussianFilterX->GetOutput() );
    g_gaussianFilterY->ReleaseDataFlagOn();
    typename g_GaussianFilterType::Pointer g_gaussianFilterZ =
      g_GaussianFilterType::New();
    g_gaussianFilterZ->SetNormalizeAcrossScale(true);
    g_gaussianFilterZ->SetSigma(currentGaussianSigma);
    g_gaussianFilterZ->SetDirection(2);
    g_gaussianFilterZ->SetInput( g_gaussianFilterY->GetOutput() );
    g_gaussianFilterZ->ReleaseDataFlagOn();
    typename itk::RescaleIntensityImageFilter<InputImageType,
                                              SmootherImageType>::Pointer
    RealToProbMapCast =
      itk::RescaleIntensityImageFilter<InputImageType, SmootherImageType>::New();
    RealToProbMapCast->SetOutputMinimum(0);
    RealToProbMapCast->SetOutputMaximum(HUNDRED_PERCENT_VALUE);
    RealToProbMapCast->SetInput( g_gaussianFilterZ->GetOutput() );
    RealToProbMapCast->Update();
    return RealToProbMapCast->GetOutput();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception in gaussian smoothing." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }
}

/**
  * \ingroup Util
  * \short This function determines the input vector size
  *        requirement based on the iris and gradient size.
  * \param NumberOfProbabilityMaps Number of Probability Maps
  * \param NumberOfImageTypes Number of Image Types (T1, T2, ..)
  * \param GradientProfileSize Desired Gradient Size
  * \param IrisSize Desired Iris Size
  * \return Vector Size Required
  */
extern unsigned int InputVectorSizeRequirement(const unsigned int NumberOfProbabilityMaps,
                                               const unsigned int NumberOfImageTypes,
                                               const unsigned int GradientProfileSize/*,
                                          *     const unsigned int IrisSize*/);

/**
  * \ingroup Util
  * \short
  *  This function compiles input vector information.
  * \param inputvector Pointer to Input Vector
  * \param CurrentIndex Current Index
  * \param Normalize Normalization Factor
  * \param MaxIndex Max Image Size
  * \param DeformedProbabilityMap Deformed Probability Map
  * \param ImageTypeList Image Type List (T1, T2, ..)
  * \param MapOfImages Map of Input Images
  * \param GradientProfileSize Desired Gradient Size
  * \param IrisSize Desired Iris Size
  * \return Void
  */
extern  void AddInputVector(std::vector<neural_scalar_type> & inputvector, itk::Image<itk::CovariantVector<float,
                                                                                                           3>,
                                                                                      3>::Pointer ProbMapGradient,
                            InternalImageType::Pointer RhoMapImage, InternalImageType::Pointer PhiMapImage,
                            InternalImageType::Pointer ThetaMapImage, const std::map<std::string,
                                                                                     InternalImageType
                                                                                     ::
                                                                                     Pointer> & MapOfImages,
                            const std::map<std::string,
                                           InternalImageType
                                           ::
                                           Pointer> & MapOfDeformedROIs, std::map<std::string,
                                                                                  neural_scalar_type>
                            & MapOfLocalMean, std::map<std::string,
                                                       neural_scalar_type>
                            & MapOfLocalStd, std::map<int,
                                                      std
                                                      ::string>& MapOfROIOrder,
                            const InternalImageType::IndexType & CurrentIndex,
                            const int GradientProfileSize);

extern int AddROIVectorTrain( ProbabilityMapParser * currentROI, DataSet * subjectSet, NetConfiguration & ANNXMLObject,
                              std::map<int,
                                       std
                                       ::
                                       string> MapOfROIOrder, std::map<std::string,
                                                                       InternalImageType
                                                                       ::
                                                                       Pointer>& MapOfImages, std::map<std::string,
                                                                                                       InternalImageType
                                                                                                       ::
                                                                                                       Pointer> &
                              MapOfDeformedROIs, std::map<std::string,
                                                          InternalImageType
                                                          ::
                                                          Pointer>& MapOfDeformedSpatialDescription,
                              const int inputVectorSize, const int outputVectorSize,
                              std::ofstream & ANNVectorStream);

extern void AddROIVectorApply( ProbabilityMapParser * currentROI, DataSet * subjectSet, NetConfiguration & ANNXMLObject,
                               std::map<int,
                                        std
                                        ::
                                        string> MapOfROIOrder, std::map<std::string,
                                                                        InternalImageType
                                                                        ::
                                                                        Pointer>& MapOfImages, std::map<std::string,
                                                                                                        InternalImageType
                                                                                                        ::
                                                                                                        Pointer> &
                               MapOfDeformedROIs, std::map<std::string,
                                                           InternalImageType
                                                           ::
                                                           Pointer>& MapOfDeformedSpatialDescription,
                               const int inputVectorSize, const int outputVectorSize,
                               const std::string ANNModelFilename);

extern int AddSubjectInputVector(DataSet * subjectSet, NetConfiguration & ANNXMLObject,
                                 const std::string registrationID, const int inputVectorSize,
                                 const int outputVectorSize, const map<int, std::string>& MapOfROIOrder,
                                 bool Apply = false);

extern void XYZToSpherical(const itk::Point<float, 3> & LocationWithOriginAtCenterOfImage, float & rho, float & phi,
                           float & theta);

/**
  * \ingroup Util
  * \short
  *  This function calculates and fills the input vector with gradient signal
  *intensity information.
  * \param fi Input Vector Iterator
  * \param MapOfImages Map of Input Images
  * \param MapOfImageInterpolators
  * \param DeformedProbMap -- Get the gradient from the
  *        probability map rather than the image since this is
  *        likely more stable.
  * \param CurrentIndex Current Index
  * \param ProfileExtent How big the profile "plug" should below
  * \param IrisSize Desired Iris Size
  * \return Void
  */
extern void FillGradProfile(std::vector<float>::iterator & fi, const std::map<std::string,
                                                                              InternalImageType::Pointer>& MapOfImages,
                            std::map<std::string,
                                     ImageLinearInterpolatorType
                                     ::
                                     Pointer>
                            MapOfImageInterpolators,
                            //  const InternalImageType::Pointer DeformedProbMap,
                            const InternalImageType::IndexType & CurrentIndex,
                            const int ProfileExtent/*,
                           * const int IrisSize*/);

/**
  * \ingroup Util
  * \short
  *  This function allocates memory for a new double image.
  * \param PointerToOutputImage Output Image Pointer (double)
  * \param PreInitializedImage Input Image Pointer (double)
  * \return Void
  */
extern void CreateNewImageFromTemplate(InternalImageType::Pointer & PointerToOutputImage,
                                       const InternalImageType::Pointer & PreInitializedImage);

/**
  * \ingroup Util
  * \short
  *  This function allocates memory for a new float image.
  * \param PointerToOutputImage Output Image Pointer (float)
  * \param PreInitializedImage Input Image Pointer (float)
  * \return Void
  */
extern void CreateNewFloatImageFromTemplate(InputImageType::Pointer & PointerToOutputImage,
                                            const InternalImageType::Pointer & PreInitializedImage);

/**
  * \ingroup Util
  * \short
  *   This function computes the bounding box of non-zero data.
  * \param image Input Image Pointer
  * \param min Minimum Indeces
  * \param max Maximum Indeces
  * \return Void
  */

extern void DefineBoundingBox(const InternalImageType::Pointer image, itk::Index<3> & min, itk::Index<3> & max);

/**
  * \ingroup Util
  * \short
  *  This function applies a landmark thin plate spline registration.
  * \param InputImageFilename Moving Image Filename
  * \param InputAtlasFilename Landmark Filename for Moving Image
  * \param TemplateAtlasFilename Landmark Filename for Fixed Image
  * \param TemplateImageFilename Fixed Image Filename
  * \param ResultingImageFilename Deformed Image Filename
  * \param MetaData Image Meta Data
  * \return Void
  */
extern void LandmarkRegistraterImages(const std::string & InputImageFilename, const std::string & InputAtlasFilename,
                                      const std::string & TemplateAtlasFilename,
                                      const std::string & TemplateImageFilename,
                                      const std::string & ResultingImageFilename, const int MetaData);

extern void PreLandmarkRegisterImage(std::string & MovingImageFilename, const std::string & PreRegisteredFilename,
                                     const std::string & LandmarkImage, const std::string & LandmarkFixedImage,
                                     const std::string & FixedImageFilename);

extern TDeformationField::Pointer RegisterLandmarksToDeformationField(const std::string & InputImageFilename,
                                                                      const std::string & InputAtlasFilename,
                                                                      const std::string & TemplateAtlasFilename,
                                                                      const std::string & TemplateImageFilename);

extern int CreateMITransformFile(const std::string & MovingImageFilename, const std::string & FixedImageFilename,
                                 const std::string & MovingImageToFixedImageRegistrationFilename,
                                 const std::string & ParamName);

extern int CreateTransformFile(const std::string & MovingImageFilename, const std::string & FixedImageFilename,
                               const std::string & OutputRegName, const std::string & FixedBinaryImageFilename,
                               const std::string & MovingBinaryImageFilename, const int roiAutoDilateSize = 1,
                               bool verbose = true);

extern int GenerateRegistrations(NetConfiguration & prob, bool reverse, bool apply, const unsigned int numThreads);

template <class WarperImageType>
typename WarperImageType::Pointer ImageWarper(
  const std::string & RegistrationFilename,
  const std::string & ImageName,
  typename WarperImageType::Pointer ReferenceImage
  )
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  const bool useTransform = ( RegistrationFilename.find(".mat") != std::string::npos );

  typename WarperImageType::Pointer PrincipalOperandImage;   // One name for the
                                                             // image to be
                                                             // warped.
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
    {
    typedef typename itk::ImageFileReader<WarperImageType> ReaderType;
    typename ReaderType::Pointer imageReader = ReaderType::New();

    imageReader->SetFileName(ImageName);
    imageReader->Update();
    PrincipalOperandImage = imageReader->GetOutput();
    }

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  typedef float                                        VectorComponentType;
  typedef typename itk::Vector<VectorComponentType, 3> VectorPixelType;
  typedef typename itk::Image<VectorPixelType,  3>     DeformationFieldType;

  // An empty SmartPointer constructor sets up someImage.IsNull() to represent
  // a
  // not-supplied state:
  typename DeformationFieldType::Pointer DeformationField;
  // typename WarperImageType::Pointer ReferenceImage;
  // if there is no *mat file.
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  if( !useTransform )    // that is, it's a warp by deformation field:
    {
    typedef typename itk::ImageFileReader<DeformationFieldType> DefFieldReaderType;
    typename DefFieldReaderType::Pointer fieldImageReader = DefFieldReaderType::New();
    fieldImageReader->SetFileName(RegistrationFilename);
    fieldImageReader->Update();
    DeformationField = fieldImageReader->GetOutput();

    // Resample deformation filed to reference image size Check that
    // deformation
    // field and reference image have same dimensions.

    // and---  ReferenceImage.IsNull() represents the delayed default
    }

  // An empty SmartPointer constructor sets up someTransform.IsNull() to
  // represent a not-supplied state:
  typename GenericTransformType::Pointer genericTransform;
  if( useTransform )    // IF there EXIST *mat file.
    {
    std::cout << "!!!!!!!!!!!! CAUTION !!!!!!!!!!!!!!!!!!!" << std::endl
              << "* Mat file exists!" << std::endl
              << "!!!!!!!!!!!! CAUTION !!!!!!!!!!!!!!!!!!!" << std::endl;
    genericTransform = itk::ReadTransformFromDisk(RegistrationFilename);
    }
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  const double defaultValue = 0;
  const typename std::string interpolationMode = "Linear";
  const typename std::string pixelType = "short";

  typename WarperImageType::Pointer TransformedImage =
    GenericTransformImage<WarperImageType, WarperImageType, DeformationFieldType>(
      PrincipalOperandImage,
      ReferenceImage,
      DeformationField,
      genericTransform,
      defaultValue,
      interpolationMode,
      pixelType == "binary");

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  return TransformedImage;
}

template <class ImageType>
typename ImageType::Pointer
ReadImage(const std::string & filename )
{
  std::cout << __FILE__ << "::" << __LINE__ << std::endl
            << filename << std::endl;

  typedef typename ImageType::Pointer ImagePointer;

  ImagePointer inputImage = itkUtil::ReadImage<ImageType>(filename);
  inputImage =
    itkUtil::ScaleAndCast<InternalImageType, InternalImageType>(inputImage,
                                                                0,
                                                                HUNDRED_PERCENT_VALUE);
  if( inputImage.IsNull() )
    {
    std::string error("Can't open image file ");
    error += filename;
    throw itk::ExceptionObject( error.c_str() );
    }
  return inputImage;
}

template <class ImageType>
typename ImageType::Pointer
ReadMedianFilteredImage(const std::string & filename,
                        const typename ImageType::SizeType & radius)

{
  typedef typename itk::MedianImageFilter<ImageType, ImageType> MIFilterType;
  typedef typename MIFilterType::Pointer                        MIFilterPointer;
  typedef typename ImageType::Pointer                           ImagePointer;
  ImagePointer inputImage = itkUtil::ReadImage<ImageType>(filename);
  inputImage =
    itkUtil::ScaleAndCast<InputImageType, InternalImageType>(inputImage,
                                                             0,
                                                             HUNDRED_PERCENT_VALUE);
  /*std::cout << "100%::: " << HUNDRED_PERCENT_VALUE << std::endl
            << filename
            << " ** Origin ** :"
            << inputImage->GetOrigin()
            << std::endl;
  */
  if( inputImage.IsNull() )
    {
    std::string error("Can't open image file ");
    error += filename;
    throw itk::ExceptionObject( error.c_str() );
    }
  MIFilterPointer filter = MIFilterType::New();
  filter->SetRadius(radius);
  filter->SetInput(inputImage);
  filter->Update();
  return filter->GetOutput();
}

/** VerifyNonZeroImage ensures that the image is not just all
  *  zeros.  The function causes the program to exit if an all
  *  zero image is encountered.
  *  \param image -- The Image to check
  *   */
template <class ImageType>
void VerifyNonZeroImage(const typename ImageType::Pointer image,
                        std::string ImageIDName)
{
  typename itk::MinimumMaximumImageCalculator<ImageType>::Pointer minmaxcalc =
    itk::MinimumMaximumImageCalculator<ImageType>::New();
  minmaxcalc->SetImage(image);
  minmaxcalc->Compute();
//  itk::Index< 3 > min, max;
//  std::cout
//  << "Minimum Value is: "
//  << static_cast< unsigned int >( minmaxcalc->GetMinimum() )
//  << " Maximum Value is: "
//  << static_cast< unsigned int >( minmaxcalc->GetMaximum() )
//  << std::endl;
//  std::cout << "Its min index position is : "
//            << minmaxcalc->GetIndexOfMinimum() << std::endl;
//  std::cout << "Its max index position is : "
//            << minmaxcalc->GetIndexOfMaximum() << std::endl;
//  min[0] = max[0];
//  min[1] = max[1];
//  min[2] = max[2];
  if( minmaxcalc->GetMaximum() == 0.0 )
    {
    std::cout << "ERROR:  Maximum deformed probability for " << ImageIDName
              << " is 0.  This is not valid." << std::endl;
    exit(-1);
    }
}

#endif
