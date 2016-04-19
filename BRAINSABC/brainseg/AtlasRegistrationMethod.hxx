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
#include "itkVersorRigid3DTransform.h"

// MI registration module
#include "AtlasRegistrationMethod.h"

#include "vnl/vnl_math.h"

#include "Log.h"

#include <fstream>
#include <sstream>
#include <iomanip>

#include "itkBRAINSROIAutoImageFilter.h"

itk::Transform<double, 3, 3>::Pointer MakeRigidIdentity(void)
{
  // Also append identity matrix for each image
  itk::VersorRigid3DTransform<double>::Pointer rigidIdentity = itk::VersorRigid3DTransform<double>::New();

  rigidIdentity->SetIdentity();
  GenericTransformType::Pointer genericTransform = rigidIdentity.GetPointer();
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
  m_SaveState(""),
  m_RestoreState(ITK_NULLPTR),
  m_DebugLevel(0)
{
  m_InputImageTissueRegion = ITK_NULLPTR;
  m_InputSpatialObjectTissueRegion = ITK_NULLPTR;
  m_AtlasToSubjectInitialTransform = ITK_NULLPTR;
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
::SetAtlasOriginalImageList(MapOfFloatImageVectors & NewAtlasList)
{
  m_AtlasOriginalImageList = NewAtlasList;
  m_AtlasToSubjectTransform = MakeRigidIdentity();
  m_DoneRegistration = false;
  m_RegistrationUpdateNeeded = true;
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetIntraSubjectOriginalImageList(MapOfFloatImageVectors &NewIntraSubjectOriginalImageList)
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

  int i = 0;
  for(auto mapOfModalImageListsIt = this->m_IntraSubjectOriginalImageList.begin();
      mapOfModalImageListsIt != this->m_IntraSubjectOriginalImageList.end();
      ++mapOfModalImageListsIt)
    {
    FloatImageVector::iterator currModeImageListIt = mapOfModalImageListsIt->second.begin();
    FloatImageVector::iterator intraImIt = this->m_IntraSubjectOriginalImageList[mapOfModalImageListsIt->first].begin();
    StringVector::iterator isNamesIt = this->m_IntraSubjectTransformFileNames[mapOfModalImageListsIt->first].begin();
    this->m_IntraSubjectTransforms[mapOfModalImageListsIt->first].clear(); //Ensure that pushing onto clean list
    while(currModeImageListIt != mapOfModalImageListsIt->second.end() )
      {
      if( itksys::SystemTools::FileExists( (*isNamesIt).c_str() ) )
        {
        try
          {
          muLogMacro(<< "Reading transform from file: "
                     << (*isNamesIt).c_str() << "." << std::endl);
          m_IntraSubjectTransforms[mapOfModalImageListsIt->first].push_back(itk::ReadTransformFromDisk((*isNamesIt).c_str()));
          }
        catch( ... )
          {
          muLogMacro(<< "Failed to read transform file caused exception." << (*isNamesIt) <<  std::endl );
          itkExceptionMacro(<< "Failed to read transform file " <<  (*isNamesIt));
          }
        }
      else if( m_ImageLinearTransformChoice == "Identity" )
        {
        muLogMacro(<< "Registering (Identity) image to key image." << std::endl);
        m_IntraSubjectTransforms[mapOfModalImageListsIt->first].push_back(MakeRigidIdentity());
        }
      else if ( (*intraImIt).GetPointer() == this->m_KeySubjectImage.GetPointer() )
        {
        muLogMacro(<< "Key image registered to itself with Identity transform." << std::endl);
        m_IntraSubjectTransforms[mapOfModalImageListsIt->first].push_back(MakeRigidIdentity());
        }
      else // when m_ImageLinearTransformChoice == "Rigid"
        {
        typedef itk::BRAINSFitHelper HelperType;
        HelperType::Pointer intraSubjectRegistrationHelper = HelperType::New();
        intraSubjectRegistrationHelper->SetSamplingPercentage(0.05); //Sample 5% of image
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
        intraSubjectRegistrationHelper->SetFixedVolume(this->GetModifiableKeySubjectImage() );
        // TODO: Find way to turn on histogram equalization for same mode images
        const int dilateSize = 15;
        const int closingSize = 15;
        intraSubjectRegistrationHelper->SetMovingVolume((*intraImIt).GetPointer());
        muLogMacro( << "Generating MovingImage Mask (Intrasubject  " << i << ")" <<  std::endl );
        typedef itk::BRAINSROIAutoImageFilter<InternalImageType, itk::Image<unsigned char, 3> > ROIAutoType;
        typename ROIAutoType::Pointer  ROIFilter = ROIAutoType::New();
        ROIFilter->SetInput((*intraImIt));
        ROIFilter->SetClosingSize(closingSize);
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

        // Delayed until first use, but only create it once.
        if( m_InputImageTissueRegion.IsNull() || m_InputSpatialObjectTissueRegion.IsNull() )
          {
          muLogMacro( << "Generating FixedImage Mask (Intrasubject)" <<  std::endl );
          typedef itk::BRAINSROIAutoImageFilter<InternalImageType, itk::Image<unsigned char, 3> > LocalROIAutoType;
          ROIFilter = LocalROIAutoType::New();
          ROIFilter->SetInput(this->GetModifiableKeySubjectImage());
          ROIFilter->SetClosingSize(closingSize);
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

        muLogMacro(<< "Registering (Rigid) image " << i << " to first image." << std::endl);
        // For better registration, several linear registration methods are run,
        // but at the end, rigid component is extracted from output linear transform.
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
        //
        // intraSubjectRegistrationHelper->SetBackgroundFillValue(backgroundFillValue);
        // NOT VALID When using initializeTransformMode
        //
        const std::string initializeTransformMode("useCenterOfHeadAlign");
        intraSubjectRegistrationHelper->SetInitializeTransformMode(initializeTransformMode);
        intraSubjectRegistrationHelper->SetMaskInferiorCutOffFromCenter(65.0); //
        //
        // maskInferiorCutOffFromCenter);
        intraSubjectRegistrationHelper->SetCurrentGenericTransform(ITK_NULLPTR);
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
        itk::VersorRigid3DTransform<double>::Pointer versorRigid = itk::ComputeRigidTransformFromGeneric(
          intraSubjectRegistrationHelper->GetCurrentGenericTransform()->GetNthTransform(0).GetPointer() );
        GenericTransformType::Pointer p = versorRigid.GetPointer();
        this->m_IntraSubjectTransforms[mapOfModalImageListsIt->first].push_back(p);
        // Write out intermodal matricies
        muLogMacro(<< "Writing " << (*isNamesIt) << "." << std::endl);
        itk::WriteTransformToDisk<double, float>(p, (*isNamesIt));
        }
      ++currModeImageListIt;
      ++isNamesIt;
      ++intraImIt;
      }
    i++;
    }
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::AverageIntraSubjectRegisteredImages()
{
  muLogMacro(<< "Warp intra subject images within one modality to the first image of that modality channel..." << std::endl);

  // Sanity Checks
  if(  m_IntraSubjectTransforms.size() != m_IntraSubjectOriginalImageList.size() )
    {
    muLogMacro( << "ERROR:  size of input image list does not match to size of intra subject transform list." << std::endl );
    muLogMacro( << "m_IntraSubjectTransforms.size(): " << m_IntraSubjectTransforms.size() << std::endl );
    muLogMacro( << "m_IntraSubjectOriginalImageList.size(): " << m_IntraSubjectOriginalImageList.size() << std::endl );
    itkExceptionMacro(<< "ERROR:  size of input image list does not match to size of intra subject transform list. "
                      << "m_IntraSubjectTransforms.size() = " << m_IntraSubjectTransforms.size() <<   std::endl
                      << "m_IntraSubjectOriginalImageList.size() = " << m_IntraSubjectOriginalImageList.size() );
    }

  this->m_RegisteredIntraSubjectImagesList.clear(); //Ensure that pushing onto clean list
  // Use resampleInPlace to transform all images to the physical space of the first key image.
  // Then, use linear interpolation to resample all within modality images to a same voxel space.
  this->m_RegisteredIntraSubjectImagesList =
    ResampleInPlaceImageList("Linear",
                             this->m_IntraSubjectOriginalImageList,
                             this->m_IntraSubjectTransforms);

  muLogMacro(<< "Average co-registered Intra subject images" << std::endl);
  this->m_ModalityAveragedOfIntraSubjectImages.clear(); //Ensure that pushing onto clean list
  this->m_ModalityAveragedOfIntraSubjectImages.resize(0); //Ensure that pushing onto clean list

  for(MapOfFloatImageVectors::iterator mapOfRegisteredModalImageListsIt = this->m_RegisteredIntraSubjectImagesList.begin();
      mapOfRegisteredModalImageListsIt != this->m_RegisteredIntraSubjectImagesList.end();
      ++mapOfRegisteredModalImageListsIt)
    {
    // If number of image of a modality (e.g. T1, T2, etc) is greater than one; then, the average image between them is computed
    const int numbOfImagesPerModality = mapOfRegisteredModalImageListsIt->second.size();
    if( numbOfImagesPerModality > 1 )
      {
      //
      this->m_ModalityAveragedOfIntraSubjectImages.push_back(
        AverageImageList<InternalImageType>(mapOfRegisteredModalImageListsIt->second)
      );
      }
    else if( numbOfImagesPerModality == 1 )
      {
      FloatImageVector::iterator intraImIt = this->m_RegisteredIntraSubjectImagesList[mapOfRegisteredModalImageListsIt->first].begin(); // each intra subject image
      this->m_ModalityAveragedOfIntraSubjectImages.push_back( (*intraImIt).GetPointer() );
      }
    else
      {
      std::cout << "NO images for modality: " << mapOfRegisteredModalImageListsIt->first << std::endl;
      }
    }
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::RegisterAtlasToSubjectImages()
{
  // Sanity Checks
  // currently we have atlases only for T1 and T2 modalities. However, we should still be able to
  // use other available input modality images (e.g. FLAIR, IDWI, PD, etc) for the segmentation process,
  // even if they will not be used in atlas to subject registration.
  if(  m_AtlasOriginalImageList.size() != m_IntraSubjectOriginalImageList.size() )
    {
    muLogMacro( << "* WARNING *:  atlas and template image list sizes do not match. " <<   std::endl );
    muLogMacro( << "There are atlas images for " << m_AtlasOriginalImageList.size() << " modalities;" <<   std::endl );
    muLogMacro( << "However, input images belong to " << m_IntraSubjectOriginalImageList.size() << " modality channels." <<   std::endl );
    muLogMacro( << "Only the first " << m_AtlasOriginalImageList.size()
                << " modality channels will be used in atlas to subject registration." <<   std::endl );
    }

  muLogMacro(<< "\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl);
  muLogMacro(<< "Register atlas to subject images..." << std::endl);
  muLogMacro(<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" << std::endl);

  /*****  Shortcut if the registration has been done previously. ******/
  // If this final transform filename exists,
  // it will be just read in and will be used directly
  // without doing the registration.
  if( itksys::SystemTools::FileExists( this->m_AtlasToSubjectTransformFileName.c_str() ) )
    {
    try
      {
      std::cout << "****************************************" << std::endl;
      muLogMacro(<< "Reading cached Atlas to subject transform: "
                 << this->m_AtlasToSubjectTransformFileName << "." << std::endl);
      m_AtlasToSubjectTransform = itk::ReadTransformFromDisk(this->m_AtlasToSubjectTransformFileName);
      // Note:  No need to write this transform to disk
      std::cout << "****************************************" << std::endl;
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
  else /***** We actaully need to run a registration ******/
    {
    typedef itk::BRAINSFitHelper HelperType;
    HelperType::Pointer atlasToSubjectRegistrationHelper = HelperType::New();
      { // Set common parameters
      atlasToSubjectRegistrationHelper->SetSamplingPercentage(0.05); //Sample 5% of image
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
    bool runSyNFull = true; // This flag runs SyN registration in mulitple resolution levels.
      {
      std::cout << "****************************************" << std::endl;
      if( m_AtlasLinearTransformChoice == "SyN" && this->m_RestoreState.IsNotNull() )
        {
        std::cout << "SyN registration is restored from state, and no atlasToSubjectInitialTransform is needed."
                  << std::endl;
        this->m_AtlasToSubjectInitialTransform = ITK_NULLPTR;
        runSyNFull = false; // When SyN is restored from state, we need to run that only in full resolution level.
        }
      else if( this->m_AtlasToSubjectInitialTransform.IsNotNull() )
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

    muLogMacro(<< "Registering first atlas image to first subject image." << std::endl);
    // Initialize the outputTransform with the initializer before starting the loop.
    this->m_AtlasToSubjectTransform = this->m_AtlasToSubjectInitialTransform;
    if( this->m_AtlasToSubjectTransform.IsNotNull() )
      {
      typedef itk::CompositeTransform<double, 3>                   CompositeTransformType;
      CompositeTransformType::Pointer atlasToSubjectCompositeTransform =
         dynamic_cast<CompositeTransformType *>( m_AtlasToSubjectTransform.GetPointer() );
      if( atlasToSubjectCompositeTransform.IsNull() )
        {
        atlasToSubjectCompositeTransform = CompositeTransformType::New();
        atlasToSubjectCompositeTransform->AddTransform( m_AtlasToSubjectTransform );
        }
      atlasToSubjectRegistrationHelper->SetCurrentGenericTransform( atlasToSubjectCompositeTransform );
      }
    // Register all atlas images to first image
    // Set the fixed and moving image
    atlasToSubjectRegistrationHelper->SetFixedVolume(this->m_ModalityAveragedOfIntraSubjectImages[0]); // by AverageIntraSubjectRegisteredImages function
    atlasToSubjectRegistrationHelper->SetMovingVolume(this->GetFirstAtlasOriginalImage());
    InternalImagePointer SecondImagePointer = this->GetSecondModalityAtlasOriginalImage("T2");
    //if ( SecondImagePointer.IsNull() ) //HACK: This is not a great solution, and it requires the atlas to have a PD image
    //{
    //    muLogMacro( << "Multimodal Registration will be run using a PD image." <<   std::endl );
    //    SecondImagePointer = this->GetSecondModalityAtlasOriginalImage("PD");
    //}
    if( this->m_ModalityAveragedOfIntraSubjectImages.size() > 1  && SecondImagePointer.IsNotNull() )
        {
        muLogMacro( << "Multimodal registration will be run using the first two modalities." <<   std::endl );
        muLogMacro( << "Number of modalities is: " << this->m_ModalityAveragedOfIntraSubjectImages.size() <<  std::endl );
        //std::cout<<this->GetSecondModalityAtlasOriginalImage("T2")<<std::endl;
        atlasToSubjectRegistrationHelper->SetFixedVolume2(this->m_ModalityAveragedOfIntraSubjectImages[1]); // by AverageIntraSubjectRegisteredImages function
        atlasToSubjectRegistrationHelper->SetMovingVolume2(SecondImagePointer);
        }
    else
        {
        std::cout<< "Multimodal Registration will NOT be run." <<   std::endl;
        muLogMacro( << "Multimodal Registration will NOT be run." <<   std::endl );
        }
    muLogMacro( << "Generating MovingImage Mask (Atlas 0)" <<   std::endl );
    const int dilateSize = 10;
    const int closingSize = 15;
    typedef itk::BRAINSROIAutoImageFilter<InternalImageType, itk::Image<unsigned char, 3> > LocalROIAutoType;
    typename LocalROIAutoType::Pointer  ROIFilter = LocalROIAutoType::New();
    ROIFilter->SetInput(this->GetFirstAtlasOriginalImage());
    ROIFilter->SetClosingSize(closingSize);
    ROIFilter->SetDilateSize(dilateSize);
    ROIFilter->Update();
    atlasToSubjectRegistrationHelper->SetMovingBinaryVolume(ROIFilter->GetSpatialObjectROI() );
    if( this->m_DebugLevel > 7 )
      {
      ByteImageType::Pointer movingMaskImage = ROIFilter->GetOutput();
      typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
      ByteWriterType::Pointer writer = ByteWriterType::New();
      writer->UseCompressionOn();

      std::ostringstream oss;
      oss << this->m_OutputDebugDir << "AtlasToSubjectRegistration_MovingMask_0.nii.gz" << std::ends;
      std::string fn = oss.str();

      writer->SetInput( movingMaskImage );
      writer->SetFileName(fn.c_str() );
      writer->Update();
      muLogMacro( << __FILE__ << " " << __LINE__ << " "  <<   std::endl );
      }

    muLogMacro( << "Generating FixedImage Mask (Subject)" <<   std::endl );
    typedef itk::BRAINSROIAutoImageFilter<InternalImageType, itk::Image<unsigned char, 3> > LocalROIAutoType;
    ROIFilter = LocalROIAutoType::New();
    ROIFilter->SetInput(this->GetModifiableKeySubjectImage());
    ROIFilter->SetClosingSize(closingSize);
    ROIFilter->SetDilateSize(dilateSize);
    ROIFilter->Update();
    atlasToSubjectRegistrationHelper->SetFixedBinaryVolume(ROIFilter->GetSpatialObjectROI() );
    if( this->m_DebugLevel > 7 )
      {
      ByteImageType::Pointer fixedMaskImage = ROIFilter->GetOutput();
      typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
      ByteWriterType::Pointer writer = ByteWriterType::New();
      writer->UseCompressionOn();

      std::ostringstream oss;
      oss << this->m_OutputDebugDir << "AtlasToSubjectRegistration_FixedMask_0.nii.gz";
      std::string fn = oss.str();

      writer->SetInput( fixedMaskImage );
      writer->SetFileName(fn.c_str() );
      writer->Update();
      muLogMacro( << __FILE__ << " " << __LINE__ << " "  <<   std::endl );
      }
    // ##########################################
    // ##########################################
    // Set up registration schedule
    // ##########################################
    // ##########################################
    {   // Now do optimal registration based on requested transform choice.
    if( m_AtlasLinearTransformChoice == "Affine" )
      {
      muLogMacro(
        << "Registering (Affine) " <<  "atlas(0) to template(0) image." << std::endl);
      std::vector<double>      minimumStepSize;
      std::vector<std::string> transformType;
      // Do full registration at first iteration level if InitialTransform not given
      if( atlasToSubjectInitialTransformName == "" )   // If no initial transform, then do full multi-step
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
    else if( m_AtlasLinearTransformChoice == "BSpline" )
      {
      muLogMacro(<< "Registering (BSpline) atlas to subject "<< std::endl);
      std::vector<double> minimumStepSize(5);
      minimumStepSize[0] = 0.00005;
      minimumStepSize[1] = 0.005;
      minimumStepSize[2] = 0.005;
      minimumStepSize[3] = 0.005;
      minimumStepSize[4] = 0.005;
      atlasToSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
      std::vector<std::string> transformType(5);
      transformType[0] = "Rigid";
      transformType[1] = "ScaleVersor3D";
      transformType[2] = "ScaleSkewVersor3D";
      transformType[3] = "Affine";
      transformType[4] = "BSpline";
      atlasToSubjectRegistrationHelper->SetTransformType(transformType);
      std::vector<int> splineGridSize(3);
      splineGridSize[0] = this->m_WarpGrid[0];
      splineGridSize[1] = this->m_WarpGrid[1];
      splineGridSize[2] = this->m_WarpGrid[2];
      atlasToSubjectRegistrationHelper->SetSplineGridSize(splineGridSize);
      // Setting max displace
      atlasToSubjectRegistrationHelper->SetMaxBSplineDisplacement(6.0);
      }
    else if( m_AtlasLinearTransformChoice == "SyN" )
      {
      muLogMacro(<< "Registering (SyN) atlas to template image." << std::endl);
      std::vector<double>      minimumStepSize;
      std::vector<std::string> transformType;

      atlasToSubjectRegistrationHelper->SetSaveState( m_SaveState );

      if( m_RestoreState.IsNotNull() ) // If RestoreState is defined only one stage of SyN is needed
        {
        atlasToSubjectRegistrationHelper->SetRestoreState( m_RestoreState );
        }
      else if( atlasToSubjectInitialTransformName == "" )   // If no initial transform, then do full multi-step
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
      else if( atlasToSubjectInitialTransformName == "AffineTransform" )   // If initial Transform is Affine, then
        // update the affine stage
        {
        minimumStepSize.push_back(0.0025);
        transformType.push_back("Affine");
        }
      else if( atlasToSubjectInitialTransformName != "SyN"
              && atlasToSubjectInitialTransformName != "BSplineTransform"
              && atlasToSubjectInitialTransformName != "CompositeTransform")
        {
        itkExceptionMacro( << "ERROR: Invalid atlasToSubjectInitialTransformName"
                           << " type for m_AtlasLinearTransformChoice of type SyN" );
        }
      minimumStepSize.push_back(0.0025);
      transformType.push_back("SyN");
      atlasToSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
      atlasToSubjectRegistrationHelper->SetTransformType(transformType);
      atlasToSubjectRegistrationHelper->SetSyNFull(runSyNFull);
      }
    else
      {
      itkExceptionMacro(<< "ERROR: Invalid atlasToSubjectInitialTransformName"
                        << " type for m_AtlasLinearTransformChoice of type "
                        << m_AtlasLinearTransformChoice);
      }
    }

/*
    if( this->m_DebugLevel > 9 && m_AtlasToSubjectTransform.IsNotNull() )
      {
      muLogMacro( << "PRE_ASSIGNMENT 0 " );
      // << transformType[0] << " first of " << transformType.size() << std::endl );
      muLogMacro(<< __FILE__ << " " << __LINE__ << " "
                 << m_AtlasToSubjectTransform->GetFixedParameters() <<   std::endl );
      muLogMacro(<< __FILE__ << " " << __LINE__ << " "
                 << m_AtlasToSubjectTransform->GetParameters() <<   std::endl );
      }
*/
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
/*
    if( this->m_DebugLevel > 9 )
      {
      muLogMacro( << "POST_ASSIGNMENT0  " );
      // << transformType[0] << " first of " << transformType.size() << std::endl );
      muLogMacro(<< __FILE__ << " " << __LINE__
                 << " " << m_AtlasToSubjectTransform->GetFixedParameters() <<   std::endl );
      muLogMacro(<< __FILE__ << " " << __LINE__
                 << " " << m_AtlasToSubjectTransform->GetParameters() <<   std::endl );
      }
*/
    // End generating the best initial transform from atlas to subject.
    muLogMacro(<< "Writing " << this->m_AtlasToSubjectTransformFileName << "." << std::endl);
    itk::WriteTransformToDisk<double, float>(m_AtlasToSubjectTransform, this->m_AtlasToSubjectTransformFileName);
    }
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::RegisterImages()
{
  // Intra subject MUST be done first
  // This function registers all the input subjects to the First image (m_KeySubjectImage),
  // and provides a map of registration transforms (m_IntraSubjectTransforms).
  /*
  NOTE: Intra subject image registration must only be rigid!
  */
  this->RegisterIntraSubjectImages();

  // This function warps all of intra subject images of one modality to the first image of that modality channel
  // using the intra subject registration transforms (m_IntraSubjectTransforms).
  // Then, it averages all images within one modality together (e.g. all T1s together and T2s together),
  // and put them in "this->m_ModalityAveragedOfIntraSubjectImages".
  // Note that this->m_ModalityAveragedOfIntraSubjectImages[0] and this->m_ModalityAveragedOfIntraSubjectImages[1]
  // (related to T1 and T2 modalities) will be used in a multi modal registration framework for atlas to subject registration.
  this->AverageIntraSubjectRegisteredImages();

  // Atlas to subject registration is done as a multi-modal registration using
  // this->m_ModalityAveragedOfIntraSubjectImages[0] and this->m_ModalityAveragedOfIntraSubjectImages[1]

  std::cout << "FirstKeyAveragedImage \n" <<  this->m_ModalityAveragedOfIntraSubjectImages[0] << std::endl;
  if( this->m_ModalityAveragedOfIntraSubjectImages.size() > 1 )
    {
    std::cout << "SecondKeyAveragedImage \n" <<  this->m_ModalityAveragedOfIntraSubjectImages[1] << std::endl;
    }
  else
    {
    std::cout << "NO SECOND KEY IMAGE" << std::endl;
    }
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
