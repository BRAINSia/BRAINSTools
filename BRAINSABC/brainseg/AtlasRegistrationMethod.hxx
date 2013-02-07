#ifndef __AtlasRegistrationMethod_hxx
#define __AtlasRegistrationMethod_hxx

#include "itkAffineTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

// MI registration module
#include "AtlasRegistrationMethod.h"
#include "RegistrationParameters.h"

#include "vnl/vnl_math.h"

#include "Log.h"

#include <fstream>
#include <sstream>
#include <iomanip>

#include "itkBRAINSROIAutoImageFilter.h"

GenericTransformType::Pointer MakeRigidIdentity(void)
{
  // Also append identity matrix for each image
  VersorRigid3DTransformType::Pointer rigidIdentity = VersorRigid3DTransformType::New();

  rigidIdentity->SetIdentity();
  GenericTransformType::Pointer genericTransform = NULL;
  genericTransform = rigidIdentity.GetPointer();
  return genericTransform;
}

template <class TOutputPixel, class TProbabilityPixel>
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::AtlasRegistrationMethod() : m_WarpGrid(3, 0),
  m_UseNonLinearInterpolation(true),
  m_DoneRegistration(false),
  m_RegistrationUpdateNeeded(true),
  m_AtlasLinearTransformChoice("Affine"),
  m_ImageLinearTransformChoice("Rigid"),
  m_DebugLevel(0)
{
  m_InputImageTissueRegion = NULL;
  m_InputSpatialObjectTissueRegion = NULL;
}

template <class TOutputPixel, class TProbabilityPixel>
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::~AtlasRegistrationMethod()
{
  m_IntraSubjectTransforms.clear();
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetSuffix(std::string suffix)
{
  m_Suffix = suffix;
  m_RegistrationUpdateNeeded = true;
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetAtlasOriginalImageList(std::vector<InternalImagePointer> & NewAtlasList)
{
  m_AtlasOriginalImageList = NewAtlasList;
  m_AtlasToSubjectTransform = MakeRigidIdentity();
  m_DoneRegistration = false;
  m_RegistrationUpdateNeeded = true;
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetIntraSubjectOriginalImageList(std::vector<InternalImagePointer> & NewIntraSubjectOriginalImageList)
{
  muLogMacro(<< "Set Intrasubject original image list" << std::endl);

  const unsigned int numIntraSubjectOriginalImages = NewIntraSubjectOriginalImageList.size();
  if( numIntraSubjectOriginalImages == 0 )
    {
    itkExceptionMacro(<< "No images specified");
    }

  m_IntraSubjectOriginalImageList = NewIntraSubjectOriginalImageList;
  m_AtlasToSubjectTransform = MakeRigidIdentity();
  // Clear previous transforms
  m_IntraSubjectTransforms.clear();
  m_IntraSubjectTransforms.resize(numIntraSubjectOriginalImages);
  for( unsigned int i = 0; i < numIntraSubjectOriginalImages; i++ )
    { // Also append identity matrix for each image
    GenericTransformType::Pointer transform = MakeRigidIdentity();
    m_IntraSubjectTransforms[i] = transform;
    }
  m_DoneRegistration = false;
  m_RegistrationUpdateNeeded = true;
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::Update()
{
  muLogMacro(<< "Update" << std::endl);

  if( m_RegistrationUpdateNeeded )
    {
    m_DoneRegistration = false;
    }
  if( !m_DoneRegistration )
    {
    this->RegisterImages();
    m_RegistrationUpdateNeeded = false;
    }
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::RegisterIntraSubjectImages()
{
  muLogMacro(<< "Register Intra subject images" << std::endl);
  for( unsigned int i = 1; i < m_IntraSubjectOriginalImageList.size(); i++ )
    {
    if( itksys::SystemTools::FileExists( this->m_IntraSubjectTransformFileNames[i].c_str() ) )
      {
      try
        {
        muLogMacro(<< "Reading transform from file: "
                   << this->m_IntraSubjectTransformFileNames[i] << "." << std::endl);
        m_IntraSubjectTransforms[i] = itk::ReadTransformFromDisk(this->m_IntraSubjectTransformFileNames[i]);
        }
      catch( ... )
        {
        muLogMacro(<< "Failed to read transform file caused exception." << this->m_IntraSubjectTransformFileNames[i]
                   <<  std::endl );
        itkExceptionMacro(<< "Failed to read transform file " <<  this->m_IntraSubjectTransformFileNames[i]);
        }
      }
    else if( m_ImageLinearTransformChoice == "Identity" )
      {
      muLogMacro(<< "Registering (Identity) image " << i + 1 << " to first image." << std::endl);
      m_IntraSubjectTransforms[i] = MakeRigidIdentity();
      }
    else
      {
      typedef itk::BRAINSFitHelper HelperType;
      HelperType::Pointer intraSubjectRegistrationHelper = HelperType::New();
      intraSubjectRegistrationHelper->SetNumberOfSamples(500000);
      intraSubjectRegistrationHelper->SetNumberOfHistogramBins(50);
      std::vector<int> numberOfIterations(1);
      numberOfIterations[0] = 1500;
      intraSubjectRegistrationHelper->SetNumberOfIterations(numberOfIterations);
      //
      //
      //
      // intraSubjectRegistrationHelper->SetMaximumStepLength(maximumStepSize);
      intraSubjectRegistrationHelper->SetTranslationScale(1000);
      intraSubjectRegistrationHelper->SetReproportionScale(1.0);
      intraSubjectRegistrationHelper->SetSkewScale(1.0);
      // Register each intrasubject image mode to first image
      intraSubjectRegistrationHelper->SetFixedVolume(m_IntraSubjectOriginalImageList[0]);
      // TODO: Find way to turn on histogram equalization for same mode images
      const int dilateSize = 15;
      intraSubjectRegistrationHelper->SetMovingVolume(m_IntraSubjectOriginalImageList[i]);
        {
        muLogMacro( << "Generating MovingImage Mask (Intrasubject  " << i << ")" <<  std::endl );
        typedef itk::BRAINSROIAutoImageFilter<InternalImageType, itk::Image<unsigned char, 3> > ROIAutoType;
        typename ROIAutoType::Pointer  ROIFilter = ROIAutoType::New();
        ROIFilter->SetInput(m_IntraSubjectOriginalImageList[i]);
        ROIFilter->SetDilateSize(dilateSize); // Only use a very small non-tissue
        // region outside of head during initial
        // runnings
        ROIFilter->Update();
        ByteImageType::Pointer movingMaskImage = ROIFilter->GetOutput();
        intraSubjectRegistrationHelper->SetMovingBinaryVolume(ROIFilter->GetSpatialObjectROI() );
        if( this->m_DebugLevel > 7 )
          {
          typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
          ByteWriterType::Pointer writer = ByteWriterType::New();
          writer->UseCompressionOn();

          std::ostringstream oss;
          oss << this->m_OutputDebugDir << "IntraSubject_MovingMask_" << i <<  ".nii.gz" << std::ends;
          std::string fn = oss.str();

          writer->SetInput( movingMaskImage );
          writer->SetFileName(fn.c_str() );
          writer->Update();
          muLogMacro( << __FILE__ << " " << __LINE__ << " "  <<  std::endl );
          }
        }
        {
        // Delayed until first use, but only create it once.
        if( m_InputImageTissueRegion.IsNull() || m_InputSpatialObjectTissueRegion.IsNull() )
          {
          muLogMacro( << "Generating FixedImage Mask (Intrasubject)" <<  std::endl );
          typedef itk::BRAINSROIAutoImageFilter<InternalImageType, itk::Image<unsigned char, 3> > ROIAutoType;
          typename ROIAutoType::Pointer  ROIFilter = ROIAutoType::New();
          ROIFilter->SetInput(m_IntraSubjectOriginalImageList[0]);
          ROIFilter->SetDilateSize(dilateSize); // Only use a very small non-tissue
          // region outside of head during
          // initial runnings
          ROIFilter->Update();
          m_InputImageTissueRegion = ROIFilter->GetOutput();
          m_InputSpatialObjectTissueRegion = ROIFilter->GetSpatialObjectROI();
          if( this->m_DebugLevel > 7 )
            {
            typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
            ByteWriterType::Pointer writer = ByteWriterType::New();
            writer->UseCompressionOn();

            std::ostringstream oss;
            oss << this->m_OutputDebugDir << "IntraSubject_FixedMask_" << 0 <<  ".nii.gz" << std::ends;
            std::string fn = oss.str();

            writer->SetInput( m_InputImageTissueRegion );
            writer->SetFileName(fn.c_str() );
            writer->Update();
            muLogMacro( << __FILE__ << " " << __LINE__ << " "  <<  std::endl );
            }
          }
        intraSubjectRegistrationHelper->SetFixedBinaryVolume(m_InputSpatialObjectTissueRegion);
        }
      if( m_ImageLinearTransformChoice == "Rigid" )
        {
        muLogMacro(<< "Registering (Rigid) image " << i << " to first image." << std::endl);
        std::vector<double> minimumStepSize(1);
        minimumStepSize[0] = 0.00005;
        intraSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
        std::vector<std::string> transformType(1);
        transformType[0] = "Rigid";
        intraSubjectRegistrationHelper->SetTransformType(transformType);
        }
      else if( m_ImageLinearTransformChoice == "Affine" )
        {
        muLogMacro(<< "Registering (Affine) image " << i << " to first image." << std::endl);
        std::vector<double> minimumStepSize(4);
        minimumStepSize[0] = 0.00005;
        minimumStepSize[1] = 0.005;
        minimumStepSize[2] = 0.005;
        minimumStepSize[3] = 0.005;
        intraSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
        std::vector<std::string> transformType(4);
        transformType[0] = "Rigid";
        transformType[1] = "ScaleVersor3D";
        transformType[2] = "ScaleSkewVersor3D";
        transformType[3] = "Affine";
        intraSubjectRegistrationHelper->SetTransformType(transformType);
        }
      else if( m_ImageLinearTransformChoice == "BSpline" )
        {
        muLogMacro(<< "Registering (BSpline) image " << i << " to first subject image." << std::endl);
        std::vector<double> minimumStepSize(5);
        minimumStepSize[0] = 0.00005;
        minimumStepSize[1] = 0.005;
        minimumStepSize[2] = 0.005;
        minimumStepSize[3] = 0.005;
        minimumStepSize[4] = 0.005;
        intraSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
        std::vector<std::string> transformType(5);
        transformType[0] = "Rigid";
        transformType[1] = "ScaleVersor3D";
        transformType[2] = "ScaleSkewVersor3D";
        transformType[3] = "Affine";
        transformType[4] = "BSpline";
        intraSubjectRegistrationHelper->SetTransformType(transformType);
        std::vector<int> splineGridSize(3);
        splineGridSize[0] = this->m_WarpGrid[0];
        splineGridSize[1] = this->m_WarpGrid[1];
        splineGridSize[2] = this->m_WarpGrid[2];
        intraSubjectRegistrationHelper->SetSplineGridSize(splineGridSize);
        // Setting max displace
        intraSubjectRegistrationHelper->SetMaxBSplineDisplacement(6.0);
        // intraSubjectRegistrationHelper->SetUseExplicitPDFDerivativesMode(useExplicitPDFDerivativesMode);
        // intraSubjectRegistrationHelper->SetUseCachingOfBSplineWeightsMode(useCachingOfBSplineWeightsMode);
        }
      else if( m_ImageLinearTransformChoice == "SyN" )
        {
        muLogMacro(<< "Registering (SyN) image " << i << " to first subject image." << std::endl);
        std::vector<double> minimumStepSize(5);
        minimumStepSize[0] = 0.00005;
        minimumStepSize[1] = 0.005;
        minimumStepSize[2] = 0.005;
        minimumStepSize[3] = 0.005;
        minimumStepSize[4] = 0.000005;
        intraSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
        std::vector<std::string> transformType(5);
        transformType[0] = "Rigid";
        transformType[1] = "ScaleVersor3D";
        transformType[2] = "ScaleSkewVersor3D";
        transformType[3] = "Affine";
        transformType[4] = "SyN";
        intraSubjectRegistrationHelper->SetTransformType(transformType);
        }
      //
      // intraSubjectRegistrationHelper->SetBackgroundFillValue(backgroundFillValue);
      // NOT VALID When using initializeTransformMode
      //
      // intraSubjectRegistrationHelper->SetCurrentGenericTransform(currentGenericTransform);
      //  intraSubjectRegistrationHelper->SetUseWindowedSinc(useWindowedSinc);
      const std::string initializeTransformMode("useCenterOfHeadAlign");
      intraSubjectRegistrationHelper->SetInitializeTransformMode(initializeTransformMode);
      intraSubjectRegistrationHelper->SetMaskInferiorCutOffFromCenter(65.0); //
      //
      // maskInferiorCutOffFromCenter);
      m_IntraSubjectTransforms[i] = NULL;
      intraSubjectRegistrationHelper->SetCurrentGenericTransform(m_IntraSubjectTransforms[i]);
      if( this->m_DebugLevel > 9 )
        {
        static unsigned int IntraSubjectRegistration = 0;
        std::stringstream   ss;
        ss << std::setw(3) << std::setfill('0') << IntraSubjectRegistration;
        intraSubjectRegistrationHelper->PrintCommandLine(true, std::string("IntraSubjectRegistration") + ss.str() );
        muLogMacro( << __FILE__ << " " << __LINE__ << " "  <<  std::endl );
        IntraSubjectRegistration++;
        }
      intraSubjectRegistrationHelper->Update();
      const unsigned int actualIterations = intraSubjectRegistrationHelper->GetActualNumberOfIterations();
      muLogMacro( << "Registration tool " << actualIterations << " iterations." << std::endl );
      m_IntraSubjectTransforms[i] = intraSubjectRegistrationHelper->GetCurrentGenericTransform();
      // Write out intermodal matricies
        {
        muLogMacro(<< "Writing " << this->m_IntraSubjectTransformFileNames[i] << "." << std::endl);
        WriteTransformToDisk(m_IntraSubjectTransforms[i], this->m_IntraSubjectTransformFileNames[i]);
        }
      }
    }
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::RegisterAtlasToSubjectImages()
{
  // Sanity Checks
  if(  (  m_AtlasOriginalImageList.size() != m_IntraSubjectTransforms.size() ) ||
       ( m_AtlasOriginalImageList.size() != m_IntraSubjectOriginalImageList.size() ) )
    {
    muLogMacro( << "ERROR:  atlas and template image list sizes do not match. " <<   std::endl );
    muLogMacro( << "m_AtlasOriginalImageList.size() = " << m_AtlasOriginalImageList.size() <<   std::endl );
    muLogMacro( << "m_IntraSubjectTransforms.size() = " << m_IntraSubjectTransforms.size() <<   std::endl );
    muLogMacro(<< "m_IntraSubjectOriginalImageList.size() = "
               << m_IntraSubjectOriginalImageList.size() <<   std::endl );
    itkExceptionMacro(<< "ERROR:  atlas and template image list sizes do not match. "
                      << "m_AtlasOriginalImageList.size() = " << m_AtlasOriginalImageList.size() <<   std::endl
                      << "m_IntraSubjectTransforms.size() = " << m_IntraSubjectTransforms.size() <<   std::endl
                      << "m_IntraSubjectOriginalImageList.size() = " << m_IntraSubjectOriginalImageList.size() );
    }

  // TODO:  Need to make this register all the atlas filenames to all the
  //       reference images.  Should probably do it in reverse order.
  muLogMacro(<< "Register atlas to subject images" << std::endl);
  if( itksys::SystemTools::FileExists( this->m_AtlasToSubjectTransformFileName.c_str() ) ) // Shortcut if the
                                                                                           // registration has been done
                                                                                           // previously.
    {
    try
      {
      muLogMacro(<< "Reading cached Atlas to subject transform: "
                 << this->m_AtlasToSubjectTransformFileName << "." << std::endl);
      m_AtlasToSubjectTransform = itk::ReadTransformFromDisk(this->m_AtlasToSubjectTransformFileName);
      // Note:  No need to write this transform to disk
      }
    catch( ... )
      {
      muLogMacro(<< "Failed to read transform file caused exception."
                 << this->m_AtlasToSubjectTransformFileName <<  std::endl );
      itkExceptionMacro(<< "Failed to read transform file " <<  this->m_AtlasToSubjectTransformFileName);
      }
    return;
    }
  else if( m_AtlasLinearTransformChoice == "Identity" )
    {
    muLogMacro(<< "Registering (Identity) atlas to first image." << std::endl);
    m_AtlasToSubjectTransform = MakeRigidIdentity();
    // Note:  No need to write this transform to disk
    return;
    }
  else // We actaully need to run a registration
    {
    typedef itk::BRAINSFitHelper HelperType;
    HelperType::Pointer atlasToSubjectRegistrationHelper = HelperType::New();
      { // Set common parameters
      atlasToSubjectRegistrationHelper->SetNumberOfSamples(500000);
      atlasToSubjectRegistrationHelper->SetNumberOfHistogramBins(50);
      std::vector<int> numberOfIterations(1);
      numberOfIterations[0] = 1500;
      atlasToSubjectRegistrationHelper->SetNumberOfIterations(numberOfIterations);
      // atlasToSubjectRegistrationHelper->SetMaximumStepLength(maximumStepSize);
      atlasToSubjectRegistrationHelper->SetTranslationScale(1000);
      atlasToSubjectRegistrationHelper->SetReproportionScale(1.0);
      atlasToSubjectRegistrationHelper->SetSkewScale(1.0);
      }
    // Deal with creating an initial transform.
    std::string atlasToSubjectInitialTransformName = "";
      {
      std::cout << "****************************************" << std::endl;
      if( this->m_AtlasToSubjectInitialTransform.IsNotNull() )
        {
        atlasToSubjectInitialTransformName = this->m_AtlasToSubjectInitialTransform->GetNameOfClass();
        std::cout << "atlasToSubjectInitialTransformName = " << atlasToSubjectInitialTransformName << std::endl;
        const std::string initializeTransformMode("Off");
        atlasToSubjectRegistrationHelper->SetInitializeTransformMode(initializeTransformMode);
        }
      else
        {
        std::cout << "No atlasToSubjectInitialTransformName specified on command line, using CenterOfHeadAlign."
                  << std::endl;
        const std::string initializeTransformMode("useCenterOfHeadAlign");
        atlasToSubjectRegistrationHelper->SetInitializeTransformMode(initializeTransformMode);
        atlasToSubjectRegistrationHelper->SetMaskInferiorCutOffFromCenter(65.0);
        }
      std::cout << "****************************************" << std::endl;
      //
      // atlasToSubjectRegistrationHelper->SetBackgroundFillValue(backgroundFillValue);
      // atlasToSubjectRegistrationHelper->SetUseWindowedSinc(useWindowedSinc);
      // TODO:  Need to make external/internal variable inside
      // Update that changes a variable to set the initialization mode "useCenterOfHeadAlign" works well for brain
      // images, but
      // will break algorithm for many other segmentation types.
      }

    muLogMacro(<< "Registering first atlas images to first subject image." << std::endl);
    // Initialize the outputTransform with the initializer before starting the loop.
    this->m_AtlasToSubjectTransform = this->m_AtlasToSubjectInitialTransform;
    atlasToSubjectRegistrationHelper->SetCurrentGenericTransform(m_AtlasToSubjectTransform);
    const unsigned int atlasReferenceImageIndex = 0;
      {
        {   // Register all atlas images to first image
          { // Set the fixed and moving image
          atlasToSubjectRegistrationHelper->SetFixedVolume(m_IntraSubjectOriginalImageList[atlasReferenceImageIndex]);
          atlasToSubjectRegistrationHelper->SetMovingVolume(m_AtlasOriginalImageList[atlasReferenceImageIndex]);
          }
          {
          muLogMacro( << "Generating MovingImage Mask (Atlas " << atlasReferenceImageIndex << ")" <<   std::endl );
          typedef itk::BRAINSROIAutoImageFilter<InternalImageType, itk::Image<unsigned char, 3> > ROIAutoType;
          typename ROIAutoType::Pointer  ROIFilter = ROIAutoType::New();
          ROIFilter->SetInput(m_AtlasOriginalImageList[atlasReferenceImageIndex]);
          ROIFilter->SetDilateSize(1); // Only use a very small non-tissue
          // region outside of head during
          // initial runnings
          ROIFilter->Update();
          atlasToSubjectRegistrationHelper->SetMovingBinaryVolume(ROIFilter->GetSpatialObjectROI() );
          if( this->m_DebugLevel > 7 )
            {
            ByteImageType::Pointer movingMaskImage = ROIFilter->GetOutput();
            typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
            ByteWriterType::Pointer writer = ByteWriterType::New();
            writer->UseCompressionOn();

            std::ostringstream oss;
            oss << this->m_OutputDebugDir << "Atlas_MovingMask_"
                << atlasReferenceImageIndex <<  ".nii.gz" << std::ends;
            std::string fn = oss.str();

            writer->SetInput( movingMaskImage );
            writer->SetFileName(fn.c_str() );
            writer->Update();
            muLogMacro( << __FILE__ << " " << __LINE__ << " "  <<   std::endl );
            }
          }
          {
          muLogMacro( << "Generating FixedImage Mask (Atlas)" <<   std::endl );
          typedef itk::BRAINSROIAutoImageFilter<InternalImageType, itk::Image<unsigned char, 3> > ROIAutoType;
          typename ROIAutoType::Pointer  ROIFilter = ROIAutoType::New();
          ROIFilter->SetInput(m_IntraSubjectOriginalImageList[atlasReferenceImageIndex]);
          ROIFilter->SetDilateSize(1); // Only use a very small non-tissue
          // region outside of head during
          // initial runnings
          ROIFilter->Update();
          atlasToSubjectRegistrationHelper->SetFixedBinaryVolume(ROIFilter->GetSpatialObjectROI() );
          if( this->m_DebugLevel > 7 )
            {
            ByteImageType::Pointer fixedMaskImage = ROIFilter->GetOutput();
            typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
            ByteWriterType::Pointer writer = ByteWriterType::New();
            writer->UseCompressionOn();

            std::ostringstream oss;
            oss << this->m_OutputDebugDir << "SubjectForAtlas_FixedMask_"
                << atlasReferenceImageIndex <<  ".nii.gz";
            std::string fn = oss.str();

            writer->SetInput( fixedMaskImage );
            writer->SetFileName(fn.c_str() );
            writer->Update();
            muLogMacro( << __FILE__ << " " << __LINE__ << " "  <<   std::endl );
            }
          }
        // ##########################################
        // ##########################################
        // Set up registration schedule
        // ##########################################
        // ##########################################
          { // Now do optimal registration based on requested transform choice.
          if( m_AtlasLinearTransformChoice == "Affine" )
            {
            muLogMacro(
              << "Registering (Affine) " <<  "atlas(" << atlasReferenceImageIndex
              << ") to template("
              << atlasReferenceImageIndex << ") image." << std::endl);
            std::vector<double>      minimumStepSize;
            std::vector<std::string> transformType;
            // Do full registration at first iteration level if InitialTransform not given
            if( atlasToSubjectInitialTransformName == "" ) // If no initial transform, then do full multi-step
                                                           // registration.
              {
              minimumStepSize.push_back(0.0025);
              transformType.push_back("Rigid");
              minimumStepSize.push_back(0.0025);
              transformType.push_back("ScaleVersor3D");
              minimumStepSize.push_back(0.0025);
              transformType.push_back("ScaleSkewVersor3D");
              }
            else if( atlasToSubjectInitialTransformName != "AffineTransform" )
              {
              itkExceptionMacro( << "ERROR: Invalid atlasToSubjectInitialTransformName"
                                 << " type for m_AtlasLinearTransformChoice of type Affine" );
              }
            minimumStepSize.push_back(0.0025);
            transformType.push_back("Affine");
            atlasToSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
            atlasToSubjectRegistrationHelper->SetTransformType(transformType);
            }
          else if( m_AtlasLinearTransformChoice == "SyN" )
            {
            muLogMacro(
              << "Registering (SyN) " << "atlas(" << atlasReferenceImageIndex
              << ") to template("
              << atlasReferenceImageIndex << ") image." << std::endl);
            std::vector<double>      minimumStepSize;
            std::vector<std::string> transformType;
            if( atlasToSubjectInitialTransformName == "" ) // If no initial transform, then do full multi-step
                                                           // registration.
              {
              minimumStepSize.push_back(0.0025);
              transformType.push_back("Rigid");
              minimumStepSize.push_back(0.0025);
              transformType.push_back("ScaleVersor3D");
              minimumStepSize.push_back(0.0025);
              transformType.push_back("ScaleSkewVersor3D");
              minimumStepSize.push_back(0.0025);
              transformType.push_back("Affine");
              }
            else if( atlasToSubjectInitialTransformName == "AffineTransform" ) // If initial Transform is Affine, then
                                                                               // update the affine stage
              {
              minimumStepSize.push_back(0.0025);
              transformType.push_back("Affine");
              }
            else if( atlasToSubjectInitialTransformName != "SyN" )
              {
              itkExceptionMacro( << "ERROR: Invalid atlasToSubjectInitialTransformName"
                                 << " type for m_AtlasLinearTransformChoice of type SyN" );
              }
            minimumStepSize.push_back(0.0025);
            transformType.push_back("SyN");
            atlasToSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
            atlasToSubjectRegistrationHelper->SetTransformType(transformType);
            }
          else
            {
            itkExceptionMacro(<< "ERROR: Invalid atlasToSubjectInitialTransformName"
                              << " type for m_AtlasLinearTransformChoice of type "
                              << m_AtlasLinearTransformChoice);
            }
          }

        if( this->m_DebugLevel > 9 && m_AtlasToSubjectTransform.IsNotNull() )
          {
          muLogMacro( << "PRE_ASSIGNMENT" <<  atlasReferenceImageIndex << "  " );
          // << transformType[0] << " first of " << transformType.size() << std::endl );
          muLogMacro(<< __FILE__ << " " << __LINE__ << " "
                     << m_AtlasToSubjectTransform->GetFixedParameters() <<   std::endl );
          muLogMacro(<< __FILE__ << " " << __LINE__ << " "
                     << m_AtlasToSubjectTransform->GetParameters() <<   std::endl );
          }
        if( this->m_DebugLevel > 9 )
          {
          static unsigned int originalAtlasToSubject = 0;
          std::stringstream   ss;
          ss << std::setw(3) << std::setfill('0') << originalAtlasToSubject;
          atlasToSubjectRegistrationHelper->PrintCommandLine(true, std::string("AtlasToSubject") + ss.str() );
          muLogMacro( << __FILE__ << " " << __LINE__ << " "  <<   std::endl );
          originalAtlasToSubject++;
          }
        atlasToSubjectRegistrationHelper->Update();
        const unsigned int actualIterations = atlasToSubjectRegistrationHelper->GetActualNumberOfIterations();
        muLogMacro( << "Registration tool " << actualIterations << " iterations." << std::endl );
        m_AtlasToSubjectTransform = atlasToSubjectRegistrationHelper->GetCurrentGenericTransform();
        if( this->m_DebugLevel > 9 )
          {
          muLogMacro( << "POST_ASSIGNMENT" <<  atlasReferenceImageIndex << "  " );
          // << transformType[0] << " first of " << transformType.size() << std::endl );
          muLogMacro(<< __FILE__ << " " << __LINE__
                     << " " << m_AtlasToSubjectTransform->GetFixedParameters() <<   std::endl );
          muLogMacro(<< __FILE__ << " " << __LINE__
                     << " " << m_AtlasToSubjectTransform->GetParameters() <<   std::endl );
          }
        }
      } // End generating the best initial transform for atlas T1 to subject T1
#if 0
    // Update the registration with all the other non-primary images.
      {
      for( unsigned int atlasIter = 1; atlasIter < m_AtlasOriginalImageList.size(); atlasIter++ )
        {
        // Save code snippets for future integration 2012-10-18
          { // Set the fixed image
          InternalImageType::Pointer currentWarpedIntraSubject = NULL;
          // NOTE: This is to save memory, so just resample the images as needed
          // to avoid keeping an extra copy of them
          if( atlasReferenceImageIndex == 0 )
            { // First subject image does not need to be resampled
            currentWarpedIntraSubject = m_IntraSubjectOriginalImageList[0];
            }
          else
            { // All other subject images must be resampled
            typedef itk::ResampleImageFilter<InternalImageType, InternalImageType> ResampleType;
            typedef ResampleType::Pointer                                          ResamplePointer;

            ResamplePointer resampler = ResampleType::New();
            resampler->SetInput(m_IntraSubjectOriginalImageList[atlasReferenceImageIndex]);
            resampler->SetTransform(m_IntraSubjectTransforms[atlasIter]);
            resampler->SetOutputParametersFromImage(m_IntraSubjectOriginalImageList[0]);
            resampler->SetDefaultPixelValue(0);
            resampler->Update();
            currentWarpedIntraSubject = resampler->GetOutput();
            }
          atlasToSubjectRegistrationHelper->SetFixedVolume(currentWarpedIntraSubject);
          }
        std::string preprocessMovingString("");
          { // Set the moving image
            // TODO:  Just turn histogram matching off at this point.
            // histogram matching often help in the registration for difficult
            // cases
            // Need to find some way to quantify how much it helps/hurts.
            // Changed so FLAIR doesn't get histogram matched to non-existing
            // atlas
            // const bool histogramMatch=true;//Setting histogram matching to true
          const bool histogramMatch = false;
          // const bool histogramMatch=(m_InputVolumeTypes[atlasIter]=="FLAIR")?false:true;
          if( histogramMatch )
            {
            typedef itk::HistogramMatchingImageFilter<InternalImageType,
                                                      InternalImageType> HistogramMatchingFilterType;
            HistogramMatchingFilterType::Pointer histogramfilter
              = HistogramMatchingFilterType::New();

            histogramfilter->SetInput( m_AtlasOriginalImageList[atlasIter]  );
            histogramfilter->SetReferenceImage( atlasToSubjectRegistrationHelper->GetFixedVolume() );

            histogramfilter->SetNumberOfHistogramLevels( 128 );
            histogramfilter->SetNumberOfMatchPoints( 2 );
            histogramfilter->ThresholdAtMeanIntensityOn();
            histogramfilter->Update();
            InternalImageType::Pointer equalizedMovingImage = histogramfilter->GetOutput();
            if( equalizedMovingImage->GetLargestPossibleRegion().GetSize() !=
                m_AtlasOriginalImageList[atlasIter]->GetLargestPossibleRegion().GetSize() )
              {
              itkExceptionMacro(<< "Histogram equalized image has wrong size." );
              }
            preprocessMovingString = "histogram equalized ";
            atlasToSubjectRegistrationHelper->SetMovingVolume(equalizedMovingImage);
            if( this->m_DebugLevel > 7 )
              {
              typedef itk::ImageFileWriter<InternalImageType> ShortWriterType;
              ShortWriterType::Pointer writer = ShortWriterType::New();
              writer->UseCompressionOn();

              std::string fn =
                this->m_OutputDebugDir + std::string("AtlasPostHistogram_") + ".nii.gz";

              writer->SetInput( equalizedMovingImage );
              writer->SetFileName(fn.c_str() );
              writer->Update();
              muLogMacro( << __FILE__ << " " << __LINE__ << " "  <<   std::endl );
              itkExceptionMacro(<< "Histogram match");
              }
            }
          else
            {
            atlasToSubjectRegistrationHelper->SetMovingVolume(m_AtlasOriginalImageList[atlasIter]);
            }
          }
        // muLogMacro( << "Refining transform with information from atlas " << atlasIter << " of " <<
        // m_AtlasOriginalImageList.size() << " to subject image " << atlasIter << "." << std::endl);
        // Code for updateing and refining the registration for all non-zero index locations.  Not implemented yet.
        }
      }
#endif
      {
      muLogMacro(<< "Writing " << this->m_AtlasToSubjectTransformFileName << "." << std::endl);
      WriteTransformToDisk(m_AtlasToSubjectTransform, this->m_AtlasToSubjectTransformFileName);
      }
    }
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::RegisterImages()
{
  // Intra subject MUST be done first
  this->RegisterIntraSubjectImages();
  // Atlas to subject MUST be done second
  this->RegisterAtlasToSubjectImages();
  m_DoneRegistration = true;
}

template <class TOutputPixel, class TProbabilityPixel>
typename AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::ProbabilityImagePointer
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::CopyProbabilityImage(
  typename AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
  ::InternalImagePointer img)
{
  muLogMacro(<< "CopyProbabilityImage" << std::endl);

  ProbabilityImagePointer outimg = ProbabilityImageType::New();
  outimg->CopyInformation(img);
  outimg->SetRegions(img->GetLargestPossibleRegion() );
  outimg->Allocate();

  typedef itk::ImageRegionIterator<InternalImageType> InternalIteratorType;
  InternalIteratorType inputIter(img, img->GetLargestPossibleRegion() );

  typedef itk::ImageRegionIterator<ProbabilityImageType>
    ProbabilityIteratorType;
  ProbabilityIteratorType outputIter(outimg,
                                     outimg->GetLargestPossibleRegion() );

  inputIter.GoToBegin();
  outputIter.GoToBegin();

  while( !inputIter.IsAtEnd() )
    {
    double p = inputIter.Get();
    if( p < 0.0 )
      {
      p = 0.0;
      }
    outputIter.Set(static_cast<ProbabilityImagePixelType>(p) );
    ++inputIter;
    ++outputIter;
    }

  return outimg;
}

#endif
