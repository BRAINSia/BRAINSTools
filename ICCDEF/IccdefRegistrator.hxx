#ifndef _IccdefRegistrator_txx
#define _IccdefRegistrator_txx

#include <sstream>
#include <string>
#include "IccdefRegistrator.h"
#include "itkCommand.h"
#include "ApplyField.h"
#include "itkStatisticsImageFilter.h"
#include "itkMetaImageIO.h"
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "vector"
#include "itkCheckerBoardImageFilter.h"
#include "itkIO.h"

#include "itkImageFileWriter.h"

#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include "itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkEuler3DTransform.h"
#include "itkResampleImageFilter.h"

#include "itkDisplacementFieldJacobianDeterminantFilter.h"
#include "itksys/SystemTools.hxx"

namespace itk
{
/*This function writes the displacement fields of the Deformation.*/
template <
  class TRealImage,
  class TOutputImage,
  class TFieldValue>
void IccdefRegistrator<TRealImage, TOutputImage,
                       TFieldValue>::WriteDisplacementComponents(TDisplacementField * df, std::string prefix)
{
//  typedef itk::Image<float, 3> OutputImageType;
  m_DefaultPixelValue = NumericTraits<PixelType>::One;

  // we use the vector index selection filter to break the deformation field
  // into x,y,z components.
  typedef itk::Image<FieldValueType,
                     3>                  ComponentImageType;
  typedef itk::VectorIndexSelectionCastImageFilter<TDisplacementField,
                                                   ComponentImageType> ComponentFilterType;

  std::string CurrentComponentFilename;
  try
    {
    char ext[3][14] = { "/xdisp.nii.gz", "/ydisp.nii.gz", "/zdisp.nii.gz"};

    typename ComponentFilterType::Pointer myComponentFilter
      = ComponentFilterType::New();
    myComponentFilter->SetInput(df);
    for( unsigned int extiter = 0; extiter < 3; extiter++ )
      {
      CurrentComponentFilename = prefix + ext[extiter];
      if( this->GetOutDebug() )
        {
        std::cout << "Writing Transform Image: "
                  << CurrentComponentFilename << std::endl;
        }

      myComponentFilter->SetIndex(extiter);

      typename ComponentImageType::Pointer DisplacementComponentImagePtr
        = myComponentFilter->GetOutput();

      itkUtil::WriteImage<ComponentImageType>(DisplacementComponentImagePtr,
                                              CurrentComponentFilename);
#if 0
      typedef itk::ImageFileWriter<ComponentImageType> FileWriterType;
      FileWriterType::Pointer DisplacementImageWriter = FileWriterType::New();
      DisplacementImageWriter->SetInput(DisplacementComponentImagePtr);
      DisplacementImageWriter->SetFileName( CurrentComponentFilename.c_str() );
      DisplacementImageWriter->Update();
#endif
      }
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "exception in file Displacement File Writer("
              << CurrentComponentFilename << ")" << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }
}

/*Constructor to initialize the parameters.*/
template <
  class TRealImage,
  class TOutputImage,
  class TFieldValue>
IccdefRegistrator<TRealImage, TOutputImage, TFieldValue>::IccdefRegistrator()
{
  // Images need to be set from the outside
  m_FixedImage = NULL;
  m_MovingImage = NULL;
  m_DisplacementField = NULL;
  m_BackwardDisplacementField = NULL;

  // Set up internal registrator with default components
  m_FixedImagePyramid = FixedImagePyramidType::New();
  m_FixedImagePyramid->UseShrinkImageFilterOff();
  m_MovingImagePyramid = MovingImagePyramidType::New();
  m_MovingImagePyramid->UseShrinkImageFilterOff();
  m_Registration = RegistrationType::New();

//  m_Registration->SetFixedImagePyramid (m_FixedImagePyramid);
//  m_Registration->SetMovingImagePyramid (m_MovingImagePyramid);

  m_DefaultPixelValue =  NumericTraits<typename RealImageType::PixelType>::Zero;
  // Setup an registration observer
  typedef SimpleMemberCommand<Self> CommandType;
  typename CommandType::Pointer command = CommandType::New();
  command->SetCallbackFunction(this, &Self::StartNewLevel);

  m_Tag = m_Registration->AddObserver(IterationEvent(), command);

  typedef VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<
      TDisplacementField, double> FieldInterpolatorType;

  typename FieldInterpolatorType::Pointer VectorInterpolator12
    = FieldInterpolatorType::New();

  typename FieldInterpolatorType::Pointer VectorInterpolator21
    = FieldInterpolatorType::New();

  m_Registration->GetFieldExpander12()->SetInterpolator(VectorInterpolator12);
  m_Registration->GetFieldExpander21()->SetInterpolator(VectorInterpolator21);

  // Default parameters
  m_NumberOfLevels = 1;

  m_NumberOfIterations = UnsignedIntArray(1);
  m_NumberOfIterations.Fill(10);
  m_OutputPrefix = "none";

  m_ForwardDisplacementFieldOutputName = "none";
  m_BackwardDisplacementFieldOutputName = "none";
  m_InitialFixedDisplacementFieldFilename = "none";
  m_InitialMovingDisplacementFieldFilename = "none";
  m_OutputJacobianImage = false;
  m_OutputDisplacement = false;
  m_OutputDisplacementField = false;

  m_UseHistogramMatching = false;
  m_OutDebug = false;

  m_InitialDisplacementField = NULL;
  m_ForwardDir = "";
  m_BackwardDir = "";
}

template <
  class TRealImage,
  class TOutputImage,
  class TFieldValue>
IccdefRegistrator<TRealImage, TOutputImage, TFieldValue>::~IccdefRegistrator()
{
  if( m_Tag )
    {
    m_Registration->RemoveObserver(m_Tag);
    }
}

/*Perform the registration of preprocessed images.*/
template <
  typename TRealImage,
  class TOutputImage,
  class TFieldValue>
void IccdefRegistrator<TRealImage, TOutputImage, TFieldValue>::Execute()
{
#if 0
  // Setup the image pyramids
  m_FixedImagePyramid->SetNumberOfLevels(m_NumberOfLevels);
  m_FixedImagePyramid->SetStartingShrinkFactors(
    m_FixedImageShrinkFactors.GetDataPointer() );

  m_MovingImagePyramid->SetNumberOfLevels(m_NumberOfLevels);
  m_MovingImagePyramid->
  SetStartingShrinkFactors( m_MovingImageShrinkFactors.GetDataPointer() );
#endif

  m_Registration->SetFixedImage(m_FixedImage);
  m_Registration->SetMovingImage(m_MovingImage);

  m_Registration->SetNumberOfLevels(m_NumberOfLevels);
  m_Registration->SetNumberOfIterations( m_NumberOfIterations.
                                         data_block() );

  m_ForwardDir = "forward";
  m_BackwardDir = "backward";

  if( this->m_OutputDisplacementField )
    {
    m_Registration->SetDisplacementFieldOutputNamePrefix(m_OutputPrefix);
    }

  // Setup the initial deformation field
  if( this->m_InitialFixedDisplacementFieldFilename != std::string("none")
      && this->m_InitialFixedDisplacementFieldFilename != std::string("") )
    {
    typedef   itk::ImageFileReader<TDisplacementField> FieldReaderType;
    typename FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName( m_InitialFixedDisplacementFieldFilename.c_str() );
    try
      {
      fieldReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Caught an ITK exception: " << std::endl;
      throw err;
      }

    m_Registration->SetInitialFixedDisplacementField(fieldReader->GetOutput() );
    }

  if( this->m_InitialMovingDisplacementFieldFilename != std::string("none")
      && this->m_InitialMovingDisplacementFieldFilename != std::string("") )
    {
    typedef   itk::ImageFileReader<TDisplacementField> FieldReaderType;
    typename FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName( m_InitialMovingDisplacementFieldFilename.c_str() );
    try
      {
      fieldReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Caught an ITK exception: " << std::endl;
      throw err;
      }

    m_Registration->SetInitialMovingDisplacementField(fieldReader->GetOutput() );
    }

#if 0
  if( this->m_InitialDisplacementField.IsNotNull() )
    {
//    m_Registration->SetInitialDisplacementField(this->m_InitialDisplacementField);
    m_Registration->SetInitialMovingDisplacementField(this->m_InitialDisplacementField);
    }
#endif

  // Perform the registration.
  try
    {
    m_Registration->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Caught an exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw err;
    throw err;
    }
  catch( ... )
    {
    std::
    cout << "Caught a non-ITK exception " << __FILE__ << " " << __LINE__
         << std::endl;
    }

  try
    {
    m_DisplacementField = m_Registration->GetOutput(0);
    m_BackwardDisplacementField = m_Registration->GetOutput(1);
    if( m_DisplacementField->GetDirection() != m_FixedImage->GetDirection() )
      {
      std::cout << "ERROR Directions don't match\n"
                << m_DisplacementField->GetDirection()
                << "\n"
                << m_FixedImage->GetDirection()
                << std::endl;
      exit(-1);
      }
    if( m_Tag )
      {
      m_Registration->RemoveObserver(m_Tag);
      m_Tag = 0;
      }
    m_Registration = NULL;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Caught an exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw err;
    }
  catch( ... )
    {
    std::
    cout << "Caught a non-ITK exception " << __FILE__ << " " << __LINE__
         << std::endl;
    }

#if 1
  // Write the output deformation fields if specified by the user.
  if( this->m_ForwardDisplacementFieldOutputName != std::string("none")
      && this->m_ForwardDisplacementFieldOutputName != std::string("") )
    {
    itkUtil::WriteImage<TDisplacementField>(m_DisplacementField,
                                            this->m_ForwardDisplacementFieldOutputName);
    if( this->GetOutDebug() )
      {
      std::cout << "---Forward deformation field has been written "
                << this->m_ForwardDisplacementFieldOutputName << "--" << std::endl;
      }
    }
  if( this->m_BackwardDisplacementFieldOutputName != std::string("none")
      && this->m_BackwardDisplacementFieldOutputName != std::string("") )
    {
    itkUtil::WriteImage<TDisplacementField>(m_BackwardDisplacementField,
                                            this->m_ForwardDisplacementFieldOutputName);
    if( this->GetOutDebug() )
      {
      std::cout << "---Backward deformation field has been written "
                << this->m_BackwardDisplacementFieldOutputName << "--" << std::endl;
      }
    }
#endif

  //  Write out the displacement fields specified by the user.
  if( this->GetOutputDisplacement() )
    {
    this->WriteDisplacementComponents(m_DisplacementField, m_ForwardDir);
    this->WriteDisplacementComponents(m_BackwardDisplacementField, m_BackwardDir);
    }

  if( this->GetOutputJacobianImage() )
    {
    typedef itk::DisplacementFieldJacobianDeterminantFilter<
        TDisplacementField, FieldValueType, TRealImage> JacobianFilterType;

    typename JacobianFilterType::Pointer m_ForwardJacobianFilter = JacobianFilterType::New();
    m_ForwardJacobianFilter->SetUseImageSpacing( true );
    m_ForwardJacobianFilter->SetInput(m_DisplacementField);
    m_ForwardJacobianFilter->UpdateLargestPossibleRegion();

    typename TRealImage::Pointer fjPtr = m_ForwardJacobianFilter->GetOutput();
    itkUtil::WriteImage<TRealImage>(fjPtr, this->m_ForwardDir + "/" + "Jacobian_forward.nii.gz");

    typename JacobianFilterType::Pointer m_BackwardJacobianFilter = JacobianFilterType::New();
    m_BackwardJacobianFilter->SetUseImageSpacing( true );
    m_BackwardJacobianFilter->SetInput(m_BackwardDisplacementField);
    m_BackwardJacobianFilter->UpdateLargestPossibleRegion();
    typename TRealImage::Pointer bjPtr = m_BackwardJacobianFilter->GetOutput();
    itkUtil::WriteImage<TRealImage>(bjPtr, this->m_BackwardDir + "/" + "Jacobian_backward.nii.gz");
    }

  if( this->m_OutputPrefix != std::string("none") )
    {
    /*Warp the image with the generated deformation field.*/
    typename RealImageType::Pointer DeformedMovingImagePtr(0);
    typename RealImageType::Pointer DeformedFixedImagePtr(0);

    typedef WarpImageFilter<RealImageType, RealImageType,
                            TDisplacementField> WarperType;
    typename WarperType::Pointer warper = WarperType::New();
    if( this->GetUseHistogramMatching() == true )
      {
      warper->SetInput( m_MovingImage );
      }
    else
      {
      warper->SetInput( m_UnNormalizedMovingImage );
      }

    warper->SetOutputSpacing( m_FixedImage->GetSpacing() );
    warper->SetOutputOrigin( m_FixedImage->GetOrigin() );
    warper->SetOutputDirection( m_FixedImage->GetDirection() );
    warper->SetDisplacementField( m_DisplacementField);
    warper->Update();
    DeformedMovingImagePtr = warper->GetOutput();

    typename WarperType::Pointer warperb = WarperType::New();
    warperb->SetInput(m_FixedImage);
    warperb->SetOutputSpacing( m_MovingImage->GetSpacing() );
    warperb->SetOutputOrigin( m_MovingImage->GetOrigin() );
    warperb->SetOutputDirection( m_MovingImage->GetDirection() );
    warperb->SetDisplacementField(m_BackwardDisplacementField);
    warperb->Update();
    DeformedFixedImagePtr = warperb->GetOutput();

    if( this->GetOutDebug() )
      {
      std::cout << "-----Direction of output warped image\n"
                << DeformedMovingImagePtr->GetDirection()
                << "\n-----Direction of deformation field\n"
                << this->m_DisplacementField->GetDirection() << std::endl;
      }
    /*Write the output image.*/
    if( this->m_OutputPrefix != std::string("none") )
      {
      typename OutputImageType::Pointer CastImageSptr
        = itkUtil::PreserveCast<RealImageType, OutputImageType>(
            DeformedMovingImagePtr);
      itkUtil::WriteImage<OutputImageType>(CastImageSptr,
                                           this->m_ForwardDir + "/" + "deformedMovingImage.nii.gz");

      if( this->GetOutDebug() )
        {
        std::cout << "---Deformed Moving Image has been written" << std::endl;
        }

      typename OutputImageType::Pointer CastImageSptr1
        = itkUtil::PreserveCast<RealImageType, OutputImageType>(DeformedFixedImagePtr);
      itkUtil::WriteImage<OutputImageType>(CastImageSptr1,
                                           this->m_BackwardDir + "/" + "deformedFixedImage.nii.gz");

      if( this->GetOutDebug() )
        {
        std::cout << "---Deformed Fixed Image has been written" << std::endl;
        }
      }
    }
}

// Print out the present registration level.
template <
  class TRealImage,
  class TOutputImage,
  class TFieldValue>
void IccdefRegistrator<TRealImage, TOutputImage, TFieldValue>::StartNewLevel()
{
  if( this->GetOutDebug() )
    {
    std::cout << "--- Starting level " << m_Registration->GetCurrentLevel()
              << std::endl;
    }
}
} // namespace itk
#endif
