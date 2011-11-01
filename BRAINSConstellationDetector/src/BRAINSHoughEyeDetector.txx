/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "BRAINSHoughEyeDetector.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
BRAINSHoughEyeDetector<TInputImage, TOutputImage>
::BRAINSHoughEyeDetector()
{
  /** Input Parameters */
  // Note the default parameter set is designed for the IXI database
  this->m_NumberOfSpheres     = 2;
  this->m_MinimumRadius       = 11.;
  this->m_MaximumRadius       = 13.;
  this->m_SigmaGradient       = 1.;
  this->m_Variance            = 1.;
  this->m_SphereRadiusRatio   = 1.;
  this->m_VotingRadiusRatio   = .5;
  this->m_Threshold           = 10.;
  this->m_OutputThreshold     = .8;
  this->m_GradientThreshold   = 0.;
  this->m_NbOfThreads         = 64;
  this->m_SamplingRatio       = .2;
  this->m_HoughEyeDetectorMode = 1;   // for T1-weighted image
  this->m_R1                  = 30;
  this->m_R2                  = 120;
  this->m_Theta               = 1.04719755; // 120 degrees 0.785398163;   //
                                            // spread angle = 90 deg
  this->m_ResultsDir          = ".";
  this->m_WritedebuggingImagesLevel = 0;

  /** Output parameters */
  this->m_AccumulatorImage    = TInputImage::New();
  this->m_RoIImage            = TInputImage::New();
  this->m_Ipd                 = 0;
  this->m_MaxInputPixelValue  = 0;
  this->m_MinInputPixelValue  = 0;
  this->m_OutputImage         = TOutputImage::New();
  this->m_Failure             = false;

  /** Internal parameters */
  this->m_VersorTransform     = VersorTransformType::New();
  this->m_InvVersorTransform  = VersorTransformType::New();
}

template <class TInputImage, class TOutputImage>
void
BRAINSHoughEyeDetector<TInputImage, TOutputImage>
::GenerateData()
{
  const InputImageConstPointer image( this->GetInput() );
  const InputRegionType        region = image->GetLargestPossibleRegion();
  const InputSizeType          size = region.GetSize();
    {
    /*
     * Set RoI of the input image
     */
    this->m_RoIImage->CopyInformation( image );
    this->m_RoIImage->SetRegions( image->GetLargestPossibleRegion() );
    this->m_RoIImage->Allocate();
    this->m_RoIImage->FillBuffer( 0 );

      {
      // A unit vector pointing to the anterior direction
      typename InputPointType::VectorType unitVectorAnterior;
      unitVectorAnterior[0] = 0.;
      unitVectorAnterior[1] = -1.;
      unitVectorAnterior[2] = 0.;

      if( ( this->m_R1 > 0 ) && ( this->m_R2 > this->m_R1 ) && ( this->m_Theta > 0 ) )
        {
        InputImageConstIterator It0( image, region );
        It0.GoToBegin();
        OutputImageIteratorWithIndex It1( this->m_RoIImage, region );
        It1.GoToBegin();

        while( !( It0.IsAtEnd() || It1.IsAtEnd() ) )
          {
          InputPointType currPt;
          image->TransformIndexToPhysicalPoint( It1.GetIndex(), currPt );

          // Center of head mass to current vector
          const typename InputPointType::VectorType CMtoCurrVec = currPt - this->m_CenterOfHeadMass;

          // posterior/anterior component of the vector
          const float CMtoCurrPA = CMtoCurrVec * unitVectorAnterior;

          // norm of CMtoCurrVec
          const float CMtoCurrNorm = CMtoCurrVec.GetNorm();

          // angle between current vector and unitVectorAnterior
          const float currTheta = acos( CMtoCurrPA / CMtoCurrNorm );

          if( ( CMtoCurrNorm > this->m_R1 ) &&
              ( CMtoCurrNorm < this->m_R2 ) &&
              ( currTheta < this->m_Theta ) )
            {
            It1.Set( It0.Get() );
            }

          // save max/min pixel value info
          if( It0.Get() < this->m_MinInputPixelValue )
            {
            this->m_MinInputPixelValue = It0.Get();
            }
          else if( It0.Get() > this->m_MaxInputPixelValue )
            {
            this->m_MaxInputPixelValue = It0.Get();
            }
          ++It0;
          ++It1;
          }
        }
      }

    /*
     * Sphere detection with Hough transform radial voting filter
     */
    HoughFilterPointer houghFilter = HoughFilterType::New();
    if( ( this->m_R1 > 0 ) && ( this->m_R2 > this->m_R1 ) && ( this->m_Theta > 0 ) )
      {
      houghFilter->SetInput( this->m_RoIImage );
      }
    else
      {
      std::cout << "Warning: RoI parameters are not valid"
                << "Set to ( default ) entire image" << std::endl;
      houghFilter->SetInput( image );
      }
    houghFilter->SetNumberOfSpheres( this->m_NumberOfSpheres );
    houghFilter->SetMinimumRadius( this->m_MinimumRadius );
    houghFilter->SetMaximumRadius( this->m_MaximumRadius );
    houghFilter->SetSigmaGradient( this->m_SigmaGradient );
    houghFilter->SetVariance( this->m_Variance );
    houghFilter->SetSphereRadiusRatio( this->m_SphereRadiusRatio );
    houghFilter->SetVotingRadiusRatio( this->m_VotingRadiusRatio );
    houghFilter->SetThreshold( this->m_Threshold );
    houghFilter->SetOutputThreshold( this->m_OutputThreshold );
    houghFilter->SetGradientThreshold( this->m_GradientThreshold );
    houghFilter->SetNbOfThreads( this->m_NbOfThreads );
    houghFilter->SetSamplingRatio( this->m_SamplingRatio );
    houghFilter->SetHoughEyeDetectorMode( this->m_HoughEyeDetectorMode );
    try
      {
      houghFilter->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Failed houghFilter " << std::endl;
      std::cerr << excep << std::endl;
      }
    catch( ... )
      {
      std::cout << "Failed on houghFilter exception occured" << std::endl;
      }
    this->m_AccumulatorImage = houghFilter->GetOutput();

    /*
     * Write debug image
     */
    if( this->m_WritedebuggingImagesLevel > 1 )
      {
        {
        // Write debug ROI image
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName( this->m_ResultsDir + "/HoughEyeROI.nii.gz" );
        writer->SetInput( this->m_RoIImage );
        writer->SetUseCompression( true );
        try
          {
          writer->Update();
          }
        catch( itk::ExceptionObject & excep )
          {
          std::cerr << "Cannot write the ROI image!" << std::endl;
          std::cerr << excep << std::endl;
          }
        }

        {
        // Write debug accumulator image
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName( this->m_ResultsDir + "/HoughEyeAccumulator.nii.gz" );
        writer->SetInput( this->m_AccumulatorImage );
        writer->SetUseCompression( true );
        try
          {
          writer->Update();
          }
        catch( itk::ExceptionObject & excep )
          {
          std::cerr << "Cannot write the ROI image!" << std::endl;
          std::cerr << excep << std::endl;
          }
        }
      }

    /*
     * Computing eye centers by finding the points with highest
     * accumulated PDF
     */
    // Get some basic information of the image

    const SpheresListType spheres = houghFilter->GetSpheres();
    if( spheres.size() < 2 )
      {
      std::cerr << "Error: The number of detected spheres is less than 2!" << std::endl;
      std::cerr << "The program will continue to run for generating some debug output for GUI corrector." << std::endl;
      std::cerr << "Make sure you set the debug level > 4." << std::endl;
      this->m_Failure = true;
      return;
      }

    InputIndexType indexEye1;
    SphereIterator itSpheres = spheres.begin();
    for( unsigned int i = 0; i < Dimension; ++i )
      {
      indexEye1[i] = static_cast<unsigned long int>( ( *itSpheres )->GetObjectToParentTransform()->GetOffset()[i] );
      }

    ++itSpheres;

    InputIndexType indexEye2;
    for( unsigned int i = 0; i < Dimension; ++i )
      {
      indexEye2[i] = static_cast<unsigned long int>( ( *itSpheres )->GetObjectToParentTransform()->GetOffset()[i] );
      }

    InputPointType physicalEye1;
    InputPointType physicalEye2;
      {
      image->TransformIndexToPhysicalPoint( indexEye1, physicalEye1 );
      image->TransformIndexToPhysicalPoint( indexEye2, physicalEye2 );
      }

    // We will determine the left and right of the eyes
    // Assuming that the degradation of the image does not
    // change their relative positions.
    unsigned int eye1OnLeft = 1;            // eye1 is on left = 1/right = 0
    if( physicalEye1[0] > physicalEye2[0] ) // eye1 is on the left
      {
      eye1OnLeft = 1;
      for( unsigned int i = 0; i < Dimension; ++i )
        {
        this->m_LE[i] = physicalEye1[i];
        this->m_RE[i] = physicalEye2[i];
        }
      }
    else
      {
      eye1OnLeft = 0;
      for( unsigned int i = 0; i < Dimension; ++i )
        {
        this->m_LE[i] = physicalEye2[i];   // eye2 is on the left
        this->m_RE[i] = physicalEye1[i];
        }
      }

    /*
     * Metrics Collection
     */
    // Metric 1: Adult Interpupillary Distance ( IPD )
    // Mean: 63mm
    // Most likely: 50mm - 75mm
    this->m_Ipd = 0.;
    for( unsigned int i = 0; i < Dimension; ++i )
      {
      this->m_Ipd += vnl_math_sqr( this->m_LE[i] - this->m_RE[i] );
      }
    this->m_Ipd = vcl_sqrt( this->m_Ipd );
    std::cout << "The resulted interpupilary distance is " << this->m_Ipd << " mm" << std::endl;

    if( this->m_Ipd < 40 or this->m_Ipd > 85 )
      {
      std::cerr << "WARNING: The distance is abnormal! Get ready to use a GUI corrector next." << std::endl;
      this->m_Failure = true;
      return;
      }

    /*
     * Align the image with eye centers
     */

    // Compute translation
    InputPointType physicalStartLocation;
    InputPointType physicalStopLocation;
      {
      InputIndexType startIndex = region.GetIndex();
      InputIndexType stopIndex;
      stopIndex[0] = startIndex[0] + size[0] - 1;
      stopIndex[1] = startIndex[1] + size[1] - 1;
      stopIndex[2] = startIndex[2] + size[2] - 1;
      image->TransformIndexToPhysicalPoint( startIndex, physicalStartLocation );
      image->TransformIndexToPhysicalPoint( stopIndex, physicalStopLocation );
      }

    // Compute the bound of image
    double bound[6];
    bound[0] =
      physicalStartLocation[0] < physicalStopLocation[0] ?
      physicalStartLocation[0] : physicalStopLocation[0];
    bound[1] =
      physicalStartLocation[0] >= physicalStopLocation[0] ?
      physicalStartLocation[0] : physicalStopLocation[0];
    bound[2] =
      physicalStartLocation[1] < physicalStopLocation[1] ?
      physicalStartLocation[1] : physicalStopLocation[1];
    bound[3] =
      physicalStartLocation[1] >= physicalStopLocation[1] ?
      physicalStartLocation[1] : physicalStopLocation[1];
    bound[4] =
      physicalStartLocation[2] < physicalStopLocation[2] ?
      physicalStartLocation[2] : physicalStopLocation[2];
    bound[5] =
      physicalStartLocation[2] >= physicalStopLocation[2] ?
      physicalStartLocation[2] : physicalStopLocation[2];

    // Space coordinate origin -> center of eye centers of input image
    VersorVectorType translation1;
    /** Center the image at a place close to the AC point.
     Here we use an approximation which is vertically at an IPD distance from center
     of eye centers, i.e., the center is at ( 0, -IPD, 0 ) in LPS coordinate system */
    VersorVectorType translation2;
    for( unsigned int i = 0; i < Dimension; ++i )
      {
      translation1[i] = 0.5 * ( this->m_LE[i]
                                + this->m_RE[i] );
      translation2[i] = 0;
      }
    translation2[1] = this->m_Ipd;

    // Compute rotation in radian
    // about +S-axis
    VersorVectorType rotation1;   // about +S-axis
    rotation1[0] = 0;
    rotation1[1] = 0;
    rotation1[2] = 1;

    // about +P-axis
    VersorVectorType rotation2;   // about +P-axis
    rotation2[0] = 0;
    rotation2[1] = 1;
    rotation2[2] = 0;

    // Note: this algorithm doesn't treat rotations about +L-axis
    this->m_RotAngle[0] = 0;

    // about +P-axis
    this->m_RotAngle[1] = vcl_atan( ( this->m_LE[2] - this->m_RE[2] )
                                    / ( this->m_LE[0] - this->m_RE[0] ) );
    // about +S-axis
    this->m_RotAngle[2] = -vcl_atan( ( this->m_LE[1] - this->m_RE[1] )
                                     / ( this->m_LE[0] - this->m_RE[0] ) );

    // Set affine tranformation
    this->m_VersorTransform->Translate( translation2 );
    this->m_VersorTransform->SetRotation( rotation1, -this->m_RotAngle[2] );
    this->m_VersorTransform->SetRotation( rotation2, -this->m_RotAngle[1] );
    this->m_VersorTransform->Translate( translation1 );

    // Get the inverse transform
    if( !this->m_VersorTransform->GetInverse( this->m_InvVersorTransform ) )
      {
      itkGenericExceptionMacro("Cannot get the inverse transform from Hough eye detector!");
      }

    /** The output image will have exact the same index contents
     but with modified image info so that the index-to-physical mapping
     makes the image in the physical space aligned */
    this->m_OutputImage->CopyInformation( image );
    this->m_OutputImage->SetRegions( region );

    this->m_OutputImage->SetOrigin( this->m_InvVersorTransform->GetMatrix()
                                    * image->GetOrigin() + this->m_InvVersorTransform->GetTranslation() );

    this->m_OutputImage->SetDirection( this->m_InvVersorTransform->GetMatrix()
                                       * image->GetDirection() );

    this->m_OutputImage->Allocate();

      {
      InputImageConstIterator It0( image, region );
      It0.GoToBegin();
      OutputImageIterator It1( this->m_OutputImage, region );
      It1.GoToBegin();
      while( !It0.IsAtEnd()  && !It1.IsAtEnd() )
        {
        It1.Set( It0.Get() );
        ++It0;
        ++It1;
        }
      }
    this->GraftOutput( this->m_OutputImage );
    }
}

template <class TInputImage, class TOutputImage>
void
BRAINSHoughEyeDetector<TInputImage, TOutputImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << "HoughEyeDetectorMode: " << this->m_HoughEyeDetectorMode << std::endl;
  os << "Number Of Spheres: " << this->m_NumberOfSpheres << std::endl;
  os << "Minimum Radius:  " << this->m_MinimumRadius << std::endl;
  os << "Maximum Radius: " << this->m_MaximumRadius << std::endl;
  os << "Derivative Scale : " << this->m_SigmaGradient << std::endl;
  os << "Accumulator Blur Variance: " << this->m_Variance << std::endl;
  os << "Sphere Radius Ratio: " << this->m_SphereRadiusRatio << std::endl;
  os << "Voting Radius Ratio: " << this->m_VotingRadiusRatio << std::endl;
  os << "Threshold: " << this->m_Threshold << std::endl;
  os << "Output Threshold : " << this->m_OutputThreshold << std::endl;
  os << "Gradient Threshold: " << this->m_GradientThreshold << std::endl;
  os << "NbOfThreads: " << this->m_NbOfThreads << std::endl;
  os << "Sampling Ratio: " << this->m_SamplingRatio << std::endl;
  os << "Interior Radius of RoI: " << this->m_R1 << std::endl;
  os << "Exterior Radius of RoI: " << this->m_R2 << std::endl;
  os << "Spread Angle of RoI: " << this->m_Theta << std::endl;
}
}
